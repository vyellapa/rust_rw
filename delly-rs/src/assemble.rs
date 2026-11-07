//! Port of shortpe.h assembleSplitReads.
//!
//! Read collection runs in parallel per chromosome; the per-SV read sets are then filled
//! sequentially in exactly the C++ order (they are libc++ unordered_sets whose iteration
//! order feeds the MSA).  MSA + consensus alignment run in parallel per SV.

use crate::config::Config;
use crate::cxx_hash::{hash_sr, LibcxxStringSet};
use crate::hts::*;
use crate::msa::msa;
use crate::split::{adjust_orientation, align_consensus};
use crate::tags::*;
use crate::util::*;
use rayon::prelude::*;
use std::collections::HashMap;

/// One split-read hit in collection order
struct Hit {
    svid: usize,
    seq: String,
    qual: u8,
}

fn collect_chromosome(c: &Config, valid: &ValidRegions, hdr: &Header, store: &HashMap<(i32, u64), i32>, svs: &[SVRecord], ref_index: i32) -> Vec<Hit> {
    let len = hdr.lens[ref_index as usize] as usize;
    let mut hits = vec![false; len];
    for &(pos, _) in store.keys() {
        hits[pos as usize] = true;
    }
    let mut out = Vec::new();
    for f in 0..c.files.len() {
        let mut bam = match IndexedBam::open(&c.files[f], Some(&c.genome)) {
            Some(b) => b,
            None => continue,
        };
        for iv in valid[ref_index as usize].ivs.clone() {
            bam.for_each_in(ref_index, iv.lower as i64, iv.upper as i64, |rec| {
                if rec.flag() & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
                    return true;
                }
                if (rec.qual() as u16) < c.min_map_qual || rec.tid() < 0 {
                    return true;
                }
                if !hits[rec.pos() as usize] {
                    return true;
                }
                let seed = hash_sr(rec.qname(), rec.flag());
                if let Some(&svid) = store.get(&(rec.pos() as i32, seed)) {
                    let sv = &svs[svid as usize];
                    if svid == sv.id {
                        let mut sequence = rec.seq_string();
                        let bp_point = if is_translocation(sv.svt) {
                            rec.tid() == sv.chr2
                        } else if sv.svt == 0 {
                            rec.pos() + 25 > sv.sv_start as i64
                        } else if sv.svt == 1 {
                            rec.pos() + 25 > sv.sv_end as i64
                        } else {
                            false
                        };
                        adjust_orientation(&mut sequence, bp_point, sv.svt);
                        out.push(Hit { svid: svid as usize, seq: String::from_utf8(sequence).unwrap(), qual: rec.qual() });
                    }
                }
                true
            });
        }
    }
    out
}

#[derive(Default)]
struct ReadSet {
    seqs: LibcxxStringSet,
    quals: Vec<u8>,
}

/// MSA + consensus alignment for one SV
fn assemble_one(c: &Config, hdr: &Header, seq: &RefSeq, snd: Option<&RefSeq>, mut sv: SVRecord, rs: &ReadSet) -> SVRecord {
    let reads: Vec<&[u8]> = rs.seqs.iter().map(|s| s.as_bytes()).collect();
    sv.consensus = String::from_utf8(msa(c, &reads)).unwrap();
    if !align_consensus(c, &hdr.lens, seq, snd, &mut sv) {
        sv.consensus.clear();
        sv.sr_support = 0;
        sv.sr_align_quality = 0.0;
    } else {
        let mut qv = rs.quals.clone();
        sv.mapq = qv.iter().map(|&q| q as i32).sum();
        sv.sr_map_quality = if !qv.is_empty() { median_vector(&mut qv) } else { 0 };
        sv.sr_support = rs.seqs.len() as i32;
    }
    sv
}

pub fn assemble_split_reads(c: &Config, valid: &ValidRegions, hdr: &Header, sr_store: &[HashMap<(i32, u64), i32>], svs: &mut [SVRecord]) {
    log("Split-read assembly");
    let fai = Faidx::load(&c.genome).expect("genome index");
    let nchr = hdr.n_targets();
    let chroms: Vec<i32> = (0..nchr).filter(|&r| !valid[r as usize].is_empty() && !sr_store[r as usize].is_empty()).collect();

    // Parallel read collection
    let snapshot: &[SVRecord] = svs;
    let collected: Vec<Vec<Hit>> = chroms.par_iter().map(|&r| collect_chromosome(c, valid, hdr, &sr_store[r as usize], snapshot, r)).collect();

    let mut tra_store: HashMap<usize, ReadSet> = HashMap::new();
    for (ci, &ref_index) in chroms.iter().enumerate() {
        // Sequential fill in C++ order
        let mut seq_store: HashMap<usize, ReadSet> = HashMap::new();
        for h in &collected[ci] {
            let cur = seq_store.get(&h.svid).map(|s| s.seqs.len()).unwrap_or(0);
            if (cur as u32) < c.max_read_per_sv {
                let target = if is_translocation(svs[h.svid].svt) { tra_store.entry(h.svid).or_default() } else { seq_store.entry(h.svid).or_default() };
                if target.seqs.insert(h.seq.clone()) {
                    target.quals.push(h.qual);
                }
            }
        }
        let seq = fai.fetch(&hdr.names[ref_index as usize], 0, hdr.lens[ref_index as usize] as i32).expect("reference sequence");
        let mut todo: Vec<usize> = Vec::new();
        for svid in 0..svs.len() {
            if is_translocation(svs[svid].svt) || svs[svid].chr != ref_index {
                continue;
            }
            let n = seq_store.get(&svid).map(|s| s.seqs.len()).unwrap_or(0);
            if n <= 1 {
                svs[svid].consensus.clear();
                svs[svid].sr_support = 0;
                svs[svid].sr_align_quality = 0.0;
                continue;
            }
            todo.push(svid);
        }
        let done: Vec<(usize, SVRecord)> = todo.par_iter().map(|&svid| (svid, assemble_one(c, hdr, &seq, None, svs[svid].clone(), &seq_store[&svid]))).collect();
        for (svid, sv) in done {
            svs[svid] = sv;
        }
    }

    // Translocations
    let empty = ReadSet::default();
    for ref2 in 0..nchr {
        if valid[ref2 as usize].is_empty() {
            continue;
        }
        let mut snd: Option<RefSeq> = None;
        for ref1 in (ref2 + 1)..nchr {
            if valid[ref1 as usize].is_empty() {
                continue;
            }
            let mut todo: Vec<usize> = Vec::new();
            for svid in 0..svs.len() {
                if !is_translocation(svs[svid].svt) || svs[svid].chr != ref1 || svs[svid].chr2 != ref2 {
                    continue;
                }
                let n = tra_store.get(&svid).map(|s| s.seqs.len()).unwrap_or(0);
                if n <= 1 {
                    svs[svid].consensus.clear();
                    svs[svid].sr_support = 0;
                    svs[svid].sr_align_quality = 0.0;
                    continue;
                }
                todo.push(svid);
            }
            if todo.is_empty() {
                continue;
            }
            let seq = fai.fetch(&hdr.names[ref1 as usize], 0, hdr.lens[ref1 as usize] as i32).expect("reference sequence");
            if snd.is_none() {
                snd = fai.fetch(&hdr.names[ref2 as usize], 0, hdr.lens[ref2 as usize] as i32);
            }
            let sndr = snd.as_ref();
            let done: Vec<(usize, SVRecord)> = todo.par_iter().map(|&svid| (svid, assemble_one(c, hdr, &seq, sndr, svs[svid].clone(), tra_store.get(&svid).unwrap_or(&empty)))).collect();
            for (svid, sv) in done {
                svs[svid] = sv;
            }
        }
    }
}
