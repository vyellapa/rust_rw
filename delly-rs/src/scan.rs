//! Port of shortpe.h scanPEandSR + mergeSort.
//!
//! Parallelisation: the C++ code runs one task per BAM file and walks chromosomes
//! sequentially.  Here every (file, chromosome) pair is an independent task.  The only
//! state that crosses chromosomes in the C++ code is (a) the per-read junction map, which is
//! sorted with a total order before use and is therefore order-insensitive, and (b) the
//! inter-chromosomal mate map `matetra`.  Each task records its inter-chromosomal mate
//! events in order and they are replayed sequentially in chromosome order afterwards,
//! which reproduces the C++ results exactly.

use crate::cluster::{cluster_pe, cluster_sr};
use crate::config::Config;
use crate::cxx_hash::{hash_pair, hash_pair_mate, hash_sr};
use crate::cxx_sort::{cpp_sort, CppLess};
use crate::hts::*;
use crate::junction::*;
use crate::tags::*;
use crate::util::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Inter-chromosomal mate events recorded during a chromosome scan
enum TraEvent {
    First { hv: u64, qual: u8, alen: i32 },
    Second { hv: u64, svt: usize, rec: BamAlignRecord, qual: u8 },
}

struct ChrScan {
    read_bp: ReadBp,
    bam_rec: Vec<Vec<BamAlignRecord>>,
    tra: Vec<TraEvent>,
    abnormal: u32,
}

/// Mate map value: (mapping quality, alignment length)
type MateMap = HashMap<u64, (u8, i32)>;

fn scan_chromosome(c: &Config, valid: &ValidRegions, lib: &LibraryInfo, file_c: usize, ref_index: i32) -> Option<ChrScan> {
    let mut bam = IndexedBam::open(&c.files[file_c], Some(&c.genome))?;
    let is_cram = c.files[file_c].ends_with("cram");
    if !is_cram && bam.mapped(ref_index) == 0 {
        return None;
    }
    let mut out = ChrScan { read_bp: HashMap::new(), bam_rec: vec![Vec::new(); 2 * DELLY_SVT_TRANS as usize], tra: Vec::new(), abnormal: 0 };
    let mut mate_map: MateMap = HashMap::new();
    for iv in valid[ref_index as usize].ivs.clone() {
        let mut last_aligned_pos: i32 = 0;
        let mut last_reads: HashSet<u64> = HashSet::new();
        bam.for_each_in(ref_index, iv.lower as i64, iv.upper as i64, |rec| {
            if rec.flag() & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP) != 0 {
                return true;
            }
            if (rec.qual() as u16) < c.min_map_qual || rec.tid() < 0 {
                return true;
            }
            let seed = hash_sr(rec.qname(), rec.flag());
            let mut rp: u32 = rec.pos() as u32;
            let mut sp: u32 = 0;
            for &cg in rec.cigar() {
                let op = cigar_op(cg);
                let ol = cigar_oplen(cg);
                if op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF {
                    sp = sp.wrapping_add(ol);
                    rp = rp.wrapping_add(ol);
                } else if op == BAM_CDEL {
                    if ol > c.min_ref_sep {
                        insert_junction(&mut out.read_bp, seed, rec, rp as i32, sp as i32, false);
                    }
                    rp = rp.wrapping_add(ol);
                    if ol > c.min_ref_sep {
                        insert_junction(&mut out.read_bp, seed, rec, rp as i32, sp as i32, true);
                    }
                } else if op == BAM_CINS {
                    if ol > c.min_ref_sep {
                        insert_junction(&mut out.read_bp, seed, rec, rp as i32, sp as i32, false);
                    }
                    sp = sp.wrapping_add(ol);
                    if ol > c.min_ref_sep {
                        insert_junction(&mut out.read_bp, seed, rec, rp as i32, sp as i32, true);
                    }
                } else if op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP {
                    let mut finalsp = sp as i32;
                    let mut scleft = false;
                    if sp == 0 {
                        finalsp += ol as i32;
                        scleft = true;
                    }
                    sp = sp.wrapping_add(ol);
                    if ol > c.min_clip {
                        insert_junction(&mut out.read_bp, seed, rec, rp as i32, finalsp, scleft);
                    }
                } else if op == BAM_CREF_SKIP {
                    rp = rp.wrapping_add(ol);
                }
            }

            // Paired-end clustering
            if rec.flag() & BAM_FPAIRED == 0 {
                return true;
            }
            if lib.median == 0 {
                return true;
            }
            if rec.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
                return true;
            }
            if rec.mtid() < 0 || rec.flag() & BAM_FMUNMAP != 0 {
                return true;
            }
            if valid[rec.mtid() as usize].is_empty() {
                return true;
            }
            if rec_translocation(rec) && (rec.qual() as u16) < c.min_tra_qual {
                return true;
            }
            let svt = isize_mapping_pos(rec, lib.max_isize_cutoff);
            if svt == -1 {
                return true;
            }
            if !c.want(svt) {
                return true;
            }
            if svt == 2 && (lib.max_isize_cutoff as i64) > rec.isize().abs() {
                return true;
            }
            if rec.pos() > last_aligned_pos as i64 {
                last_reads.clear();
                last_aligned_pos = rec.pos() as i32;
            }
            if first_pair_obs(rec, &last_reads) {
                last_reads.insert(seed);
                let hv = hash_pair(rec.qname(), rec.tid(), rec.pos(), rec.mtid(), rec.mpos());
                let alen = alignment_length(rec) as i32;
                if is_translocation(svt) {
                    out.tra.push(TraEvent::First { hv, qual: rec.qual(), alen });
                } else {
                    mate_map.insert(hv, (rec.qual(), alen));
                }
            } else {
                let hv = hash_pair_mate(rec.qname(), rec.tid(), rec.pos(), rec.mtid(), rec.mpos());
                let alen = alignment_length(rec);
                let mk = |pair_quality: u8, alenmate: i32| BamAlignRecord {
                    tid: rec.tid(),
                    pos: rec.pos() as i32,
                    mtid: rec.mtid(),
                    mpos: rec.mpos() as i32,
                    alen: alen as u16 as i32,
                    malen: alenmate as u16 as i32,
                    median: lib.median,
                    mad: lib.mad,
                    max_normal_isize: lib.max_normal_isize,
                    flag: rec.flag() as u32,
                    map_quality: pair_quality,
                };
                if is_translocation(svt) {
                    out.tra.push(TraEvent::Second { hv, svt: svt as usize, rec: mk(0, 0), qual: rec.qual() });
                } else {
                    match mate_map.get(&hv) {
                        Some(&(q, alenmate)) if q != 0 => {
                            mate_map.remove(&hv);
                            let pq = q.min(rec.qual());
                            out.bam_rec[svt as usize].push(mk(pq, alenmate));
                            out.abnormal += 1;
                        }
                        _ => {}
                    }
                }
            }
            true
        });
    }
    Some(out)
}

/// scanPEandSR; returns (pe svs, sr svs, srStore)
pub fn scan_pe_and_sr(c: &Config, valid: &ValidRegions, hdr: &Header, sample_lib: &mut [LibraryInfo]) -> (Vec<SVRecord>, Vec<SVRecord>, Vec<HashMap<(i32, u64), i32>>) {
    log("Paired-end and split-read scanning");
    let nf = c.files.len();
    let nchr = hdr.n_targets();
    let tasks: Vec<(usize, i32)> = (0..nf).flat_map(|f| (0..nchr).filter(|&r| !valid[r as usize].is_empty()).map(move |r| (f, r))).collect();
    let libs: Vec<LibraryInfo> = sample_lib.to_vec();
    let results: Vec<((usize, i32), Option<ChrScan>)> = tasks.par_iter().map(|&(f, r)| ((f, r), scan_chromosome(c, valid, &libs[f], f, r))).collect();

    // Per-file sequential replay in chromosome order
    let nsvt = 2 * DELLY_SVT_TRANS as usize;
    let mut per_file_sr: Vec<Vec<Vec<SRBamRecord>>> = Vec::with_capacity(nf);
    let mut per_file_pe: Vec<Vec<Vec<BamAlignRecord>>> = Vec::with_capacity(nf);
    let mut by_file: Vec<Vec<ChrScan>> = (0..nf).map(|_| Vec::new()).collect();
    for ((f, _r), res) in results.into_iter() {
        if let Some(s) = res {
            by_file[f].push(s);
        }
    }
    let file_outputs: Vec<(Vec<Vec<SRBamRecord>>, Vec<Vec<BamAlignRecord>>, u32)> = by_file
        .into_par_iter()
        .map(|scans| {
            let mut matetra: MateMap = HashMap::new();
            let mut read_bp: ReadBp = HashMap::new();
            let mut pe: Vec<Vec<BamAlignRecord>> = vec![Vec::new(); nsvt];
            let mut abnormal: u32 = 0;
            for s in scans {
                for (svt, v) in s.bam_rec.into_iter().enumerate() {
                    pe[svt].extend(v);
                }
                abnormal += s.abnormal;
                for ev in s.tra {
                    match ev {
                        TraEvent::First { hv, qual, alen } => {
                            matetra.insert(hv, (qual, alen));
                        }
                        TraEvent::Second { hv, svt, mut rec, qual } => {
                            if let Some(&(q, alenmate)) = matetra.get(&hv) {
                                if q != 0 {
                                    matetra.remove(&hv);
                                    rec.map_quality = q.min(qual);
                                    rec.malen = alenmate as u16 as i32;
                                    pe[svt].push(rec);
                                    abnormal += 1;
                                }
                            }
                        }
                    }
                }
                for (k, v) in s.read_bp {
                    read_bp.entry(k).or_default().extend(v);
                }
            }
            for v in read_bp.values_mut() {
                // Junction::operator< is a total order, so any sort gives libc++'s result
                v.sort_unstable_by(|a, b| if a.lt(b) { std::cmp::Ordering::Less } else if b.lt(a) { std::cmp::Ordering::Greater } else { std::cmp::Ordering::Equal });
            }
            let mut sr: Vec<Vec<SRBamRecord>> = vec![Vec::new(); nsvt];
            if c.want(2) {
                select_deletions(c, &read_bp, &mut sr);
            }
            if c.want(3) {
                select_duplications(c, &read_bp, &mut sr);
            }
            if c.svtset.is_empty() || c.svtset.contains(&0) || c.svtset.contains(&1) {
                select_inversions(c, &read_bp, &mut sr);
            }
            if c.want(4) {
                select_insertions(c, &read_bp, &mut sr);
            }
            if c.svtset.is_empty() || (0..4).any(|k| c.svtset.contains(&(DELLY_SVT_TRANS + k))) {
                select_translocations(c, &read_bp, &mut sr);
            }
            (sr, pe, abnormal)
        })
        .collect();
    for (f, (sr, pe, abn)) in file_outputs.into_iter().enumerate() {
        sample_lib[f].abnormal_pairs += abn;
        per_file_sr.push(sr);
        per_file_pe.push(pe);
    }
    let mut sr_br: Vec<Vec<SRBamRecord>> = vec![Vec::new(); nsvt];
    let mut bam_record: Vec<Vec<BamAlignRecord>> = vec![Vec::new(); nsvt];
    for f in 0..nf {
        for svt in 0..nsvt {
            bam_record[svt].append(&mut per_file_pe[f][svt]);
            sr_br[svt].append(&mut per_file_sr[f][svt]);
        }
    }

    // Split-read clustering
    log("Split-read clustering");
    let mut sr_svs: Vec<SVRecord> = Vec::new();
    for svt in 0..nsvt {
        if !c.want(svt as i32) || sr_br[svt].is_empty() {
            continue;
        }
        cpp_sort(&mut sr_br[svt]);
        cluster_sr(c, &mut sr_br[svt], &mut sr_svs, svt as i32);
    }

    // Paired-end clustering (independent per SV type)
    log("Paired-end clustering");
    let varisize = get_variability(sample_lib) as u32;
    let pe_parts: Vec<Vec<SVRecord>> = bam_record
        .par_iter_mut()
        .enumerate()
        .map(|(svt, recs)| {
            let mut out = Vec::new();
            if !c.want(svt as i32) || recs.is_empty() {
                return out;
            }
            cpp_sort(recs);
            cluster_pe(c, recs, &mut out, varisize, svt as i32);
            out
        })
        .collect();
    let svs: Vec<SVRecord> = pe_parts.into_iter().flatten().collect();

    // Track split-reads
    let mut sr_store: Vec<HashMap<(i32, u64), i32>> = vec![HashMap::new(); nchr as usize];
    for svt in 0..nsvt {
        for r in &sr_br[svt] {
            if r.svid != -1 && r.rstart != -1 {
                if r.rstart < hdr.lens[r.chr as usize] as i32 {
                    sr_store[r.chr as usize].entry((r.rstart, r.id)).or_insert(r.svid);
                }
                if r.chr != r.chr2 && r.rstart < hdr.lens[r.chr2 as usize] as i32 {
                    sr_store[r.chr2 as usize].entry((r.rstart, r.id)).or_insert(r.svid);
                }
            }
        }
    }
    (svs, sr_svs, sr_store)
}

/// mergeSort
pub fn merge_sort(pe: &mut Vec<SVRecord>, sr: &mut [SVRecord]) {
    cpp_sort(pe);
    cpp_sort(sr);
    for svt in 0..10 {
        for i in 0..sr.len() {
            if sr[i].svt != svt {
                continue;
            }
            if sr[i].sr_support == 0 || sr[i].sr_align_quality == 0.0 {
                continue;
            }
            let search_window = 500;
            let mut sv_exists = false;
            let key = SVRecord::key(sr[i].chr, 0.max(sr[i].sv_start - search_window), sr[i].sv_end);
            let mut k = pe.partition_point(|x| x.lt(&key));
            while k < pe.len() && (pe[k].sv_start - sr[i].sv_start).abs() < search_window {
                let o = &mut pe[k];
                k += 1;
                if o.svt != svt || o.precise {
                    continue;
                }
                if sr[i].chr != o.chr || sr[i].chr2 != o.chr2 {
                    continue;
                }
                if o.sv_start + o.ciposlow < sr[i].sv_start && sr[i].sv_start < o.sv_start + o.ciposhigh && o.sv_end + o.ciendlow < sr[i].sv_end && sr[i].sv_end < o.sv_end + o.ciendhigh {
                    sv_exists = true;
                    let s = &sr[i];
                    o.sv_start = s.sv_start;
                    o.sv_end = s.sv_end;
                    o.ciposlow = s.ciposlow;
                    o.ciposhigh = s.ciposhigh;
                    o.ciendlow = s.ciendlow;
                    o.ciendhigh = s.ciendhigh;
                    o.sr_map_quality = s.sr_map_quality;
                    o.sr_support = s.sr_support;
                    o.ins_len = s.ins_len;
                    o.hom_len = s.hom_len;
                    o.sr_align_quality = s.sr_align_quality;
                    o.precise = true;
                    o.consensus = s.consensus.clone();
                    o.cons_bp = s.cons_bp;
                    o.mapq += s.mapq;
                }
            }
            if !sv_exists {
                let prec_window = 10;
                let mut precise_dup = false;
                let ii = i as i32;
                for j in (i + 1)..sr.len() {
                    if (sr[i].sv_start - sr[j].sv_start).abs() > prec_window {
                        break;
                    }
                    if sr[i].svt != sr[j].svt || sr[i].chr != sr[j].chr || sr[i].chr2 != sr[j].chr2 {
                        continue;
                    }
                    if sr[j].sv_start + sr[j].ciposlow <= sr[i].sv_start
                        && sr[i].sv_start <= sr[j].sv_start + sr[j].ciposhigh
                        && sr[j].sv_end + sr[j].ciendlow <= sr[i].sv_end
                        && sr[i].sv_end <= sr[j].sv_end + sr[j].ciendhigh
                        && (sr[i].sr_support < sr[j].sr_support || (ii < j as i32 && sr[i].sr_support == sr[j].sr_support))
                    {
                        precise_dup = true;
                    }
                }
                for j in (0..i).rev() {
                    if (sr[i].sv_start - sr[j].sv_start).abs() > prec_window {
                        break;
                    }
                    if sr[i].svt != sr[j].svt || sr[i].chr != sr[j].chr || sr[i].chr2 != sr[j].chr2 {
                        continue;
                    }
                    if sr[j].sv_start + sr[j].ciposlow < sr[i].sv_start
                        && sr[i].sv_start < sr[j].sv_start + sr[j].ciposhigh
                        && sr[j].sv_end + sr[j].ciendlow < sr[i].sv_end
                        && sr[i].sv_end < sr[j].sv_end + sr[j].ciendhigh
                        && (sr[i].sr_support < sr[j].sr_support || (ii < j as i32 && sr[i].sr_support == sr[j].sr_support))
                    {
                        precise_dup = true;
                    }
                }
                if !precise_dup {
                    pe.push(sr[i].clone());
                    cpp_sort(pe);
                }
            }
        }
    }
}
