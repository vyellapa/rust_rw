//! Port of coverage.h: REF/ALT probe generation and SV genotyping read annotation.
//!
//! Samples are processed in parallel (all per-sample state is independent).  Within a
//! sample, ordered candidate logs are admitted in C++ batches (131072 * threads).
//! Independent SVs and read waves within an SV align in parallel. Results retain read
//! order, and caps are applied before expensive alignment as well as during admission.

use crate::config::Config;
use crate::cxx_hash::{hash_pair, hash_pair_mate, hash_string};
use crate::cxx_sort::{cpp_sort, CppLess};
use crate::edlib::{self, Mode, Task};
use crate::hts::*;
use crate::split::{adjust_orientation, cons_ref_alignment, find_split, get_sv_ref, AlignDescriptor};
use crate::svanno::annotate_sv;
use crate::tags::*;
use crate::util::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
struct SpanPoint {
    bppos: i32,
    svt: i32,
    id: u32,
    chr2: i32,
    other_bppos: i32,
}

impl CppLess for SpanPoint {
    #[inline]
    fn lt(&self, o: &Self) -> bool {
        self.bppos < o.bppos
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BpRegion {
    region_start: i32,
    region_end: i32,
    bppos: i32,
    hom_left: i32,
    hom_right: i32,
    svt: i32,
    id: u32,
    bp_point: u8,
}

impl CppLess for BpRegion {
    #[inline]
    fn lt(&self, o: &Self) -> bool {
        self.bppos < o.bppos
    }
}

#[derive(Clone, Debug, Default)]
pub struct SpanningCount {
    pub ref_: Vec<u8>,
    pub alt: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct JunctionCount {
    pub ref_: Vec<u8>,
    pub alt: Vec<u8>,
    pub hp1ref: Vec<u8>,
    pub hp1alt: Vec<u8>,
    pub hp2ref: Vec<u8>,
    pub hp2alt: Vec<u8>,
    pub ps: i32,
    pub gl: [f32; 3],
    pub locus_gq: i32,
    pub locus_gt: i8,
    pub joint: bool,
}

impl Default for JunctionCount {
    fn default() -> Self {
        JunctionCount { ref_: Vec::new(), alt: Vec::new(), hp1ref: Vec::new(), hp1alt: Vec::new(), hp2ref: Vec::new(), hp2alt: Vec::new(), ps: -1, gl: [0.0; 3], locus_gq: 0, locus_gt: -1, joint: false }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ReadCount {
    pub left_rc: i32,
    pub rc: i32,
    pub right_rc: i32,
}

/// _editDistanceHW
#[inline]
fn edit_distance_hw(c: &Config, query: &[u8], target: &[u8]) -> (f64, i32) {
    let mut score = 0.0f64;
    let mut ed = query.len() as i32;
    let k = (2.0f32 * c.flank_quality * query.len() as f32) as i32;
    let r = edlib::align(query, target, k, Mode::HW, Task::Distance);
    if r.edit_distance != -1 {
        score = ((1.0 - c.flank_quality as f64) * query.len() as f64) / (r.edit_distance + 1) as f64;
        ed = r.edit_distance;
    }
    (score, ed)
}

/// _editDistGenoQual
#[inline]
fn edit_dist_geno_qual(c: &Config, ref_ed: i32, alt_ed: i32, read_qual: u8) -> u8 {
    let delta = ref_ed - alt_ed;
    let adelta = delta.abs();
    let w = (c.flank_quality as f64 / (1.0 - c.flank_quality as f64)).log10();
    let mut ex = adelta as f64 * w;
    if ex > 4.0 {
        ex = 4.0;
    }
    let mut mq = (10.0 * (1.0 + 10f64.powf(ex)).log10()) as u32;
    if mq > c.geno_cap as u32 {
        mq = c.geno_cap as u32;
    }
    mq.min(read_qual as u32) as u8
}

fn cut_ref_start(rs: i32, re: i32, off: i32, bp_point: u32, svt: i32) -> i32 {
    let swap = if is_translocation(svt) { span_orientation(svt) == 3 } else { svt == 3 };
    if swap {
        if bp_point == 0 {
            re - off
        } else {
            rs - off
        }
    } else if bp_point != 0 {
        re - off
    } else {
        rs - off
    }
}

fn cut_ref_end(rs: i32, re: i32, off: i32, bp_point: u32, svt: i32) -> i32 {
    let swap = if is_translocation(svt) { span_orientation(svt) == 3 } else { svt == 3 };
    if swap {
        if bp_point == 0 {
            re + off
        } else {
            rs + off
        }
    } else if bp_point != 0 {
        re + off
    } else {
        rs + off
    }
}

/// std::string::substr(pos, n) with C++ int→size_t conversions
fn substr_i(s: &[u8], pos: i32, n: i32) -> Vec<u8> {
    let p = pos as usize;
    assert!(p <= s.len(), "basic_string::substr out of range");
    let cnt = if n < 0 { usize::MAX } else { n as usize };
    let end = if cnt > s.len() - p { s.len() } else { p + cnt };
    s[p..end].to_vec()
}

/// _refAnchor
fn ref_anchor(seq: &RefSeq, sv_start: i32, chr_len: i32) -> String {
    let mut p = sv_start - 1;
    if p < 0 {
        p = 0;
    }
    if p >= chr_len {
        p = chr_len - 1;
    }
    (seq.at(p as usize).to_ascii_uppercase() as char).to_string()
}

/// _addAlleles(ref, chr2, sv, svt)
fn add_alleles(r: &str, chr2: &str, sv: &SVRecord, svt: i32) -> String {
    if is_translocation(svt) {
        match span_orientation(svt) {
            0 => format!("{},{}]{}:{}]", r, r, chr2, sv.sv_end),
            1 => format!("{},[{}:{}[{}", r, chr2, sv.sv_end, r),
            2 => format!("{},{}[{}:{}[", r, r, chr2, sv.sv_end),
            3 => format!("{},]{}:{}]{}", r, chr2, sv.sv_end, r),
            _ => format!("{},<{}>", r, add_id(svt)),
        }
    } else {
        format!("{},<{}>", r, add_id(svt))
    }
}

pub struct Probes {
    pub ref_probe: [Vec<Vec<u8>>; 2],
    pub cons_probe: [Vec<Vec<u8>>; 2],
    pub bp_region: Vec<Vec<BpRegion>>,
    pub sv_on_chr: Vec<bool>,
}

/// _generateProbes
fn generate_probes(c: &Config, hdr: &Header, svs: &mut [SVRecord]) -> Probes {
    log("Generate REF and ALT probes");
    let n = svs.len();
    let nchr = hdr.n_targets() as usize;
    let mut pr = Probes { ref_probe: [vec![Vec::new(); n], vec![Vec::new(); n]], cons_probe: [vec![Vec::new(); n], vec![Vec::new(); n]], bp_region: vec![Vec::new(); nchr], sv_on_chr: vec![false; nchr] };
    let mut tra_ref: Vec<Vec<u8>> = vec![Vec::new(); n];
    let fai = Faidx::load(&c.genome).expect("genome index");
    for ref_index in 0..nchr as i32 {
        let chr_svs: Vec<usize> = (0..n).filter(|&i| svs[i].chr == ref_index || svs[i].chr2 == ref_index).collect();
        if chr_svs.is_empty() {
            continue;
        }
        pr.sv_on_chr[ref_index as usize] = true;
        let seq = fai.fetch(&hdr.names[ref_index as usize], 0, hdr.lens[ref_index as usize] as i32).expect("reference sequence");

        // Per-SV work in parallel; results applied in SV order (= single-threaded C++ order)
        struct Out {
            idx: usize,
            sv: SVRecord,
            tra_ref: Option<Vec<u8>>,
            probes: Vec<(usize, Vec<u8>, Vec<u8>, i32, BpRegion)>,
        }
        let tra_ref_ro = &tra_ref;
        let svs_ro: &[SVRecord] = svs;
        let outs: Vec<Out> = chr_svs
            .par_iter()
            .map(|&i| {
                let mut sv = svs_ro[i].clone();
                let mut out = Out { idx: i, sv: SVRecord::default(), tra_ref: None, probes: Vec::new() };
                if sv.chr == ref_index {
                    if sv.alleles.is_empty() {
                        sv.alleles = add_alleles(&ref_anchor(&seq, sv.sv_start, hdr.lens[ref_index as usize] as i32), &hdr.names[sv.chr2 as usize], &sv, sv.svt);
                    }
                    if !is_translocation(sv.svt) {
                        annotate_sv(c, &hdr.lens, &seq, &mut sv);
                    }
                }
                if sv.precise {
                    if sv.chr != sv.chr2 && sv.chr2 == ref_index {
                        let mut bp = Breakpoint::from_sv(&sv);
                        init_breakpoint(&hdr.lens, &mut bp, sv.consensus.len() as i32, sv.svt);
                        out.tra_ref = Some(get_sv_ref(c, &seq, &bp, ref_index, sv.svt));
                    }
                    if sv.chr == ref_index {
                        let mut bp = Breakpoint::from_sv(&sv);
                        if is_translocation(sv.svt) {
                            bp.part1 = tra_ref_ro[sv.id as usize].clone();
                        }
                        if sv.svt == 4 {
                            let buffer = ((sv.consensus.len() as u64).wrapping_sub(sv.ins_len as i64 as u64) / 3) as i32;
                            init_breakpoint(&hdr.lens, &mut bp, buffer.max(c.minimum_flank_size), sv.svt);
                        } else {
                            init_breakpoint(&hdr.lens, &mut bp, sv.consensus.len() as i32, sv.svt);
                        }
                        let sv_ref = get_sv_ref(c, &seq, &bp, ref_index, sv.svt);
                        let cons = sv.consensus.as_bytes();
                        if let Some(al) = cons_ref_alignment(cons, &sv_ref, sv.svt) {
                            let mut ad = AlignDescriptor::default();
                            if find_split(c, cons, &sv_ref, &al, &mut ad, sv.svt) {
                                let mfs = c.minimum_flank_size;
                                for bp_point in 0..2u32 {
                                    let (region_chr, region_start, region_end, ccs, cce, crs, cre, bppos);
                                    if bp_point == 1 {
                                        region_chr = sv.chr2;
                                        region_start = 0.max(sv.sv_end - mfs);
                                        region_end = ((sv.sv_end + mfs) as u32).min(hdr.lens[sv.chr2 as usize]) as i32;
                                        ccs = ad.c_end - ad.hom_left - mfs;
                                        cce = ad.c_end + ad.hom_right + mfs;
                                        bppos = sv.sv_end;
                                    } else {
                                        region_chr = sv.chr;
                                        region_start = 0.max(sv.sv_start - mfs);
                                        region_end = ((sv.sv_start + mfs) as u32).min(hdr.lens[sv.chr as usize]) as i32;
                                        ccs = ad.c_start - ad.hom_left - mfs;
                                        cce = ad.c_start + ad.hom_right + mfs;
                                        bppos = sv.sv_start;
                                    }
                                    crs = cut_ref_start(ad.r_start, ad.r_end, ad.hom_left + mfs, bp_point, sv.svt);
                                    cre = cut_ref_end(ad.r_start, ad.r_end, ad.hom_right + mfs, bp_point, sv.svt);
                                    let cp = substr_i(cons, ccs, cce - ccs);
                                    let rp = substr_i(&sv_ref, crs, cre - crs);
                                    let reg = BpRegion { region_start, region_end, bppos, hom_left: ad.hom_left, hom_right: ad.hom_right, svt: sv.svt, id: sv.id as u32, bp_point: bp_point as u8 };
                                    out.probes.push((bp_point as usize, cp, rp, region_chr, reg));
                                }
                            }
                        }
                    }
                }
                out.sv = sv;
                out
            })
            .collect();
        for o in outs {
            let id = o.sv.id as usize;
            if let Some(t) = o.tra_ref {
                tra_ref[id] = t;
            }
            for (bpp, cp, rp, region_chr, reg) in o.probes {
                pr.cons_probe[bpp][id] = cp;
                pr.ref_probe[bpp][id] = rp;
                pr.bp_region[region_chr as usize].push(reg);
            }
            svs[o.idx] = o.sv;
        }
    }
    for r in pr.bp_region.iter_mut() {
        cpp_sort(r);
    }
    pr
}

struct DumpRead {
    file: usize,
    header: Arc<Header>,
    qname: String,
    tid: i32,
    pos: i64,
    mtid: i32,
    mpos: i64,
    qual: u8,
}

impl DumpRead {
    fn line(&self, c: &Config, svt: i32, id: u32) -> String {
        format!("{}{:08}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tSR\n",
            add_id(svt), id, c.files[self.file], self.qname, self.header.names[self.tid as usize],
            self.pos, self.header.names[self.mtid as usize], self.mpos, self.qual)
    }
}

struct AlignJob {
    id: u32,
    bp_point: u8,
    svt: i32,
    sequence: Arc<[u8]>,
    qual: u8,
    dump: Option<Arc<DumpRead>>,
}

struct AlignResult {
    kind: u8,
    qual: u8,
    dump: Option<String>,
}

fn align_job(c: &Config, pr: &Probes, job: AlignJob) -> AlignResult {
    let id = job.id as usize;
    let mut sequence = job.sequence.to_vec();
    adjust_orientation(&mut sequence, job.bp_point != 0, job.svt);
    let (score_alt, alt_ed) = edit_distance_hw(c, &pr.cons_probe[job.bp_point as usize][id], &sequence);
    let (score_ref, ref_ed) = edit_distance_hw(c, &pr.ref_probe[job.bp_point as usize][id], &sequence);
    if score_ref > 0.7 || score_alt > 0.7 {
        let kind = if score_ref > score_alt { b'R' } else { b'A' };
        AlignResult {
            kind,
            qual: edit_dist_geno_qual(c, ref_ed, alt_ed, job.qual),
            dump: if kind == b'A' { job.dump.map(|read| read.line(c, job.svt, job.id)) } else { None },
        }
    } else {
        AlignResult { kind: b'N', qual: 0, dump: None }
    }
}

/// Ordered per-chromosome genotyping events, replayed sequentially per sample so that
/// cross-chromosome state (translocation mates, read caps, parity subsampling, value order)
/// evolves exactly as in the sequential C++ implementation.
enum Event {
    /// Junction-read candidate (one C++ AlignJob if it passes the read-cap pre-check)
    Job(AlignJob),
    /// Replaces a job after its result has been applied; emitted in original event order.
    JunctionDump(Option<String>),
    /// Normal spanning pair (reference support, subject to parity subsampling)
    SpanRef { id: u32, q: u8 },
    /// Abnormal spanning pair (variant support)
    SpanAlt { id: u32, q: u8, dump: Option<String> },
    /// First read of an inter-chromosomal pair
    TraFirst { hv: u64, clip: bool, q: u8 },
    /// Second read of an inter-chromosomal pair; `alts` are the spanning SVs it supports
    TraSecond { hv: u64, q: u8, alts: Vec<(u32, Option<String>)> },
}

struct ChromResult {
    events: Vec<Event>,
    rc: Vec<(usize, ReadCount)>,
}

#[allow(clippy::too_many_arguments)]
fn dump_line(c: &Config, file_c: usize, fhdr: &Header, svt: i32, id: u32, rec: &Record, kind: &str) -> String {
    format!(
        "{}{:08}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        add_id(svt),
        id,
        c.files[file_c],
        String::from_utf8_lossy(rec.qname()),
        fhdr.names[rec.tid() as usize],
        rec.pos(),
        fhdr.names[rec.mtid() as usize],
        rec.mpos(),
        rec.qual() as i32,
        kind
    )
}

/// Per-base yes/no mask at 1 bit per base (C++ uses one `bool` per base).
struct BitSet(Vec<u64>);

impl BitSet {
    fn new(n: usize) -> Self {
        BitSet(vec![0; n.div_ceil(64)])
    }
    #[inline]
    fn set(&mut self, i: usize) {
        self.0[i >> 6] |= 1 << (i & 63);
    }
    #[inline]
    fn get(&self, i: usize) -> bool {
        (self.0[i >> 6] >> (i & 63)) & 1 != 0
    }
}

/// Saturation limit of the C++ uint16_t coverage counters (`maxCoverage - 1`).
const COV_SATURATION: u16 = u16::MAX - 1;

/// Half-open range summed by the C++ read-count loop
/// `for (k = a; k < b && k < len; ++k)` with unsigned `k`; `None` when it is empty.
fn rc_range(a: i32, b: i32, len: u32) -> Option<(u32, u32)> {
    let (lo, hi) = (a as u32, (b as u32).min(len));
    if lo < hi {
        Some((lo, hi))
    } else {
        None
    }
}

/// Left, middle and right read-count ranges of one SV, and whether it uses base coverage
/// (small SVs) rather than fragment-midpoint coverage.
fn rc_ranges(c: &Config, sv: &SVRecord, len: u32) -> (bool, [Option<(u32, u32)>; 3]) {
    let mut small = false;
    let mut half = (sv.sv_end - sv.sv_start) / 2;
    let tra_or_ins = is_translocation(sv.svt) || sv.svt == 4;
    if tra_or_ins {
        half = 500;
        small = true;
    } else if (sv.sv_end - sv.sv_start) <= c.indelsize {
        small = true;
    }
    let left = rc_range((sv.sv_start - half).max(0), sv.sv_start, len);
    let mid = if tra_or_ins { rc_range((sv.sv_start - half).max(0), (sv.sv_start + half).min(len as i32), len) } else { rc_range(sv.sv_start, sv.sv_end, len) };
    let right = if tra_or_ins { rc_range(sv.sv_start, (sv.sv_start + half).min(len as i32), len) } else { rc_range(sv.sv_end, (sv.sv_end + half).min(len as i32), len) };
    (small, [left, mid, right])
}

/// Saturating per-base coverage, stored only inside the windows that are summed later
/// (C++ keeps a uint16_t per base of the chromosome).
struct WindowedCounts {
    /// Sorted, disjoint windows: (start, end, offset into `counts`)
    windows: Vec<(u32, u32, usize)>,
    counts: Vec<u16>,
}

impl WindowedCounts {
    fn new(mut ranges: Vec<(u32, u32)>) -> Self {
        ranges.sort_unstable();
        let mut windows: Vec<(u32, u32, usize)> = Vec::new();
        let mut total = 0usize;
        for (s, e) in ranges {
            if let Some(last) = windows.last_mut() {
                if s <= last.1 {
                    if e > last.1 {
                        total += (e - last.1) as usize;
                        last.1 = e;
                    }
                    continue;
                }
            }
            windows.push((s, e, total));
            total += (e - s) as usize;
        }
        WindowedCounts { windows, counts: vec![0; total] }
    }

    /// Count one aligned base at every position of [start, end)
    fn add(&mut self, start: u32, end: u32) {
        let mut i = self.windows.partition_point(|w| w.1 <= start);
        while i < self.windows.len() && self.windows[i].0 < end {
            let (ws, we, off) = self.windows[i];
            let (lo, hi) = (start.max(ws), end.min(we));
            for v in &mut self.counts[off + (lo - ws) as usize..off + (hi - ws) as usize] {
                if *v < COV_SATURATION {
                    *v += 1;
                }
            }
            i += 1;
        }
    }

    /// Sum over [lo, hi), which lies inside one window by construction
    fn sum(&self, lo: u32, hi: u32) -> i32 {
        let i = self.windows.partition_point(|w| w.1 <= lo);
        let (ws, _, off) = self.windows[i];
        self.counts[off + (lo - ws) as usize..off + (hi - ws) as usize].iter().fold(0i32, |s, &v| s.wrapping_add(v as i32))
    }
}

/// Saturating fragment-midpoint coverage as sorted distinct positions with prefix sums
/// (one entry per read pair instead of a uint16_t per base).
struct MidpointCounts {
    pos: Vec<u32>,
    prefix: Vec<i64>,
}

impl MidpointCounts {
    fn new(mut mids: Vec<u32>) -> Self {
        mids.sort_unstable();
        let mut pos = Vec::new();
        let mut prefix = vec![0i64];
        let mut i = 0;
        while i < mids.len() {
            let p = mids[i];
            let mut j = i;
            while j < mids.len() && mids[j] == p {
                j += 1;
            }
            pos.push(p);
            prefix.push(prefix[prefix.len() - 1] + (j - i).min(COV_SATURATION as usize) as i64);
            i = j;
        }
        MidpointCounts { pos, prefix }
    }

    /// Sum over [lo, hi); the i64 difference truncated to i32 equals the C++ wrapping int32 sum
    fn sum(&self, lo: u32, hi: u32) -> i32 {
        let a = self.pos.partition_point(|&p| p < lo);
        let b = self.pos.partition_point(|&p| p < hi);
        (self.prefix[b] - self.prefix[a]) as i32
    }
}

/// Scan one chromosome of one sample (runs in parallel)
#[allow(clippy::too_many_arguments)]
fn annotate_chrom(c: &Config, lib: &LibraryInfo, svs: &[SVRecord], pr: &Probes, hdr: &Header, file_c: usize, ref_index: i32) -> Option<ChromResult> {
    let mut bam = IndexedBam::open(&c.files[file_c], Some(&c.genome))?;
    bam.set_required_fields(SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_RNEXT | SAM_PNEXT | SAM_TLEN | SAM_SEQ | SAM_QNAME);
    let is_cram = c.files[file_c].ends_with("cram");
    if !is_cram && bam.mapped(ref_index) == 0 {
        return None;
    }
    let fhdr = Arc::new(bam.header.clone());
    let mfs = c.minimum_flank_size;
    let mut events: Vec<Event> = Vec::new();
    let mut clip: HashMap<u64, (bool, u8)> = HashMap::new();
    let tlen = fhdr.lens[ref_index as usize];
    // Read-count ranges of the SVs on this chromosome: base coverage is only kept where it is summed.
    let len0 = hdr.lens[ref_index as usize];
    let sv_ranges: Vec<(usize, bool, [Option<(u32, u32)>; 3])> = svs
        .iter()
        .filter(|sv| sv.chr == ref_index)
        .map(|sv| {
            let (small, r) = rc_ranges(c, sv, len0);
            (sv.id as usize, small, r)
        })
        .collect();
    let mut cov_bases = WindowedCounts::new(sv_ranges.iter().filter(|s| s.1).flat_map(|s| s.2.iter().flatten().copied()).collect());
    let mut mids: Vec<u32> = Vec::new();
    let regions = &pr.bp_region[ref_index as usize];
    let mut bp_occupied = BitSet::new(tlen as usize);
    for r in regions {
        for k in r.region_start..r.region_end {
            bp_occupied.set(k as usize);
        }
    }
    let mut span_point: Vec<SpanPoint> = Vec::new();
    let mut span_bp = BitSet::new(tlen as usize);
    for sv in svs {
        if sv.pe_support == 0 {
            continue;
        }
        if sv.chr == ref_index && sv.sv_start < tlen as i32 {
            span_bp.set(sv.sv_start as usize);
            span_point.push(SpanPoint { bppos: sv.sv_start, svt: sv.svt, id: sv.id as u32, chr2: sv.chr2, other_bppos: sv.sv_end });
        }
        if sv.chr2 == ref_index && sv.sv_end < tlen as i32 {
            span_bp.set(sv.sv_end as usize);
            span_point.push(SpanPoint { bppos: sv.sv_end, svt: sv.svt, id: sv.id as u32, chr2: sv.chr, other_bppos: sv.sv_start });
        }
    }
    cpp_sort(&mut span_point);

    let mut last_aligned_pos: i32 = 0;
    let mut last_reads: HashSet<u64> = HashSet::new();
    bam.for_each_in(ref_index, 0, tlen as i64, |rec| {
        if rec.flag() & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP) != 0 {
            return true;
        }
        if (rec.qual() as u16) < c.min_geno_qual {
            return true;
        }
        let pos = rec.pos();
        {
            let mut rp: u32 = 0;
            for &cg in rec.cigar() {
                let op = cigar_op(cg);
                if op == BAM_CMATCH {
                    let start = pos + rp as i64;
                    let end = (start + cigar_oplen(cg) as i64).min(tlen as i64);
                    if start < end {
                        cov_bases.add(start as u32, end as u32);
                    }
                    rp += cigar_oplen(cg);
                } else if op == BAM_CDEL || op == BAM_CREF_SKIP {
                    rp += cigar_oplen(cg);
                }
            }
        }
        let mut has_soft_clip = false;
        let mut has_clip = false;
        let mut leading_sc: i32 = 0;
        for (i, &cg) in rec.cigar().iter().enumerate() {
            let op = cigar_op(cg);
            if op == BAM_CSOFT_CLIP {
                has_clip = true;
                has_soft_clip = true;
                if i == 0 {
                    leading_sc = cigar_oplen(cg) as i32;
                }
            } else if op == BAM_CHARD_CLIP {
                has_clip = true;
            }
        }
        let lq = rec.l_qseq() as i64;
        if c.max_geno_read_count != 0 && rec.l_qseq() >= 2 * mfs {
            let rbegin: i32 = 0.max(pos as i32 - leading_sc);
            let mut bpvalid = false;
            let mut k = rbegin as i64;
            while k < pos + lq && k < tlen as i64 {
                if bp_occupied.get(k as usize) {
                    bpvalid = true;
                    break;
                }
                k += 1;
            }
            if bpvalid {
                let mut it = regions.partition_point(|b| b.bppos < rbegin);
                let mut raw: Option<Arc<[u8]>> = None;
                let mut dump_read: Option<Arc<DumpRead>> = None;
                while it < regions.len() && pos + lq >= regions[it].bppos as i64 {
                    let b = regions[it];
                    it += 1;
                    if has_soft_clip || (!has_clip && pos + (mfs + b.hom_left) as i64 <= b.bppos as i64 && pos + lq >= (b.bppos + mfs + b.hom_right) as i64) {
                        let sequence = Arc::clone(raw.get_or_insert_with(|| Arc::from(rec.seq_string())));
                        let dump = if c.has_dump_file {
                            Some(Arc::clone(dump_read.get_or_insert_with(|| Arc::new(DumpRead {
                                file: file_c, header: Arc::clone(&fhdr), qname: String::from_utf8_lossy(rec.qname()).into_owned(),
                                tid: rec.tid(), pos: rec.pos(), mtid: rec.mtid(), mpos: rec.mpos(), qual: rec.qual(),
                            }))))
                        } else { None };
                        events.push(Event::Job(AlignJob { id: b.id, bp_point: b.bp_point, svt: b.svt, sequence, qual: rec.qual(), dump }));
                    }
                }
            }
        }

        // Read-count and spanning annotation
        if rec.flag() & BAM_FPAIRED == 0 || !pr.sv_on_chr[rec.mtid() as usize] {
            return true;
        }
        if pos > last_aligned_pos as i64 {
            last_reads.clear();
            last_aligned_pos = pos as i32;
        }
        let intra = rec.tid() == rec.mtid();
        if first_pair_obs(rec, &last_reads) {
            last_reads.insert(hash_string(rec.qname()) as u64);
            let hv = hash_pair(rec.qname(), rec.tid(), rec.pos(), rec.mtid(), rec.mpos());
            if intra {
                clip.insert(hv, (has_soft_clip, rec.qual()));
            } else {
                events.push(Event::TraFirst { hv, clip: has_soft_clip, q: rec.qual() });
            }
            return true;
        }
        let hv = hash_pair_mate(rec.qname(), rec.tid(), rec.pos(), rec.mtid(), rec.mpos());
        if !intra {
            // Inter-chromosomal: the mate lookup is deferred to the replay.  Such a pair can only
            // contribute abnormal spanning support (the other branches require tid == mtid).
            let mut alts: Vec<(u32, Option<String>)> = Vec::new();
            if lib.median != 0 {
                let svt = isize_mapping_pos(rec, lib.max_isize_cutoff);
                if svt != -1 {
                    let mut pbegin = pos as i32;
                    let mut pend = (pos as i32 + lib.max_normal_isize).min(tlen as i32);
                    if rec.flag() & BAM_FREVERSE != 0 {
                        pbegin = 0.max(pos as i32 + rec.l_qseq() - lib.max_normal_isize);
                        pend = (pos as i32 + rec.l_qseq()).min(tlen as i32);
                    }
                    if (pbegin..pend).any(|i| span_bp.get(i as usize)) {
                        let mut it = span_point.partition_point(|s| s.bppos < pbegin);
                        while it < span_point.len() && pend >= span_point[it].bppos {
                            let sp = span_point[it];
                            it += 1;
                            if svt == sp.svt && rec.mtid() == sp.chr2 && (rec.mpos() as i32 - sp.other_bppos).abs() < lib.max_normal_isize {
                                let d = if c.has_dump_file { Some(dump_line(c, file_c, &fhdr, sp.svt, sp.id, rec, "PE")) } else { None };
                                alts.push((sp.id, d));
                            }
                        }
                    }
                }
            }
            events.push(Event::TraSecond { hv, q: rec.qual(), alts });
            return true;
        }
        let (pair_quality, pair_clip) = match clip.remove(&hv) {
            None => return true,
            Some((sc, q)) => (q.min(rec.qual()), sc || has_soft_clip),
        };
        if (pair_quality as u16) < c.min_geno_qual {
            return true;
        }
        let mid = pos as i32 + half_alignment_length(rec) as i32;
        if mid < tlen as i32 {
            mids.push(mid as u32);
        }
        let outer: i32 = if rec.pos() < rec.mpos() { (rec.mpos() + lq - rec.pos()) as i32 } else { (rec.pos() + lq - rec.mpos()) as i32 };
        if lib.median == 0 {
            return true;
        }
        let svtype = get_sv_type(rec);
        if !pair_clip && svtype == 2 && outer >= lib.min_normal_isize && outer <= lib.max_normal_isize {
            let spanlen = (0.8 * outer as f64) as i32;
            let pbegin = (rec.pos() as i32).min(rec.mpos() as i32);
            let st = pbegin + (outer - spanlen) / 2;
            let mut valid = false;
            let mut i = st;
            while i < st + spanlen && i < tlen as i32 {
                if span_bp.get(i as usize) {
                    valid = true;
                    break;
                }
                i += 1;
            }
            if valid {
                let mut it = span_point.partition_point(|s| s.bppos < st);
                while it < span_point.len() && st + spanlen >= span_point[it].bppos {
                    events.push(Event::SpanRef { id: span_point[it].id, q: pair_quality });
                    it += 1;
                }
            }
        }
        if svtype != 2 || outer < lib.min_normal_isize || outer > lib.max_normal_isize {
            let svt = isize_mapping_pos(rec, lib.max_isize_cutoff);
            if svt == -1 {
                return true;
            }
            let mut pbegin = pos as i32;
            let mut pend = (pos as i32 + lib.max_normal_isize).min(tlen as i32);
            if rec.flag() & BAM_FREVERSE != 0 {
                pbegin = 0.max(pos as i32 + rec.l_qseq() - lib.max_normal_isize);
                pend = (pos as i32 + rec.l_qseq()).min(tlen as i32);
            }
            if (pbegin..pend).any(|i| span_bp.get(i as usize)) {
                let mut it = span_point.partition_point(|s| s.bppos < pbegin);
                while it < span_point.len() && pend >= span_point[it].bppos {
                    let sp = span_point[it];
                    it += 1;
                    if svt == sp.svt && rec.mtid() == sp.chr2 && (rec.mpos() as i32 - sp.other_bppos).abs() < lib.max_normal_isize {
                        let d = if c.has_dump_file { Some(dump_line(c, file_c, &fhdr, sp.svt, sp.id, rec, "PE")) } else { None };
                        events.push(Event::SpanAlt { id: sp.id, q: pair_quality, dump: d });
                    }
                }
            }
        }
        true
    });

    // Assign fragment and base counts to SVs on this chromosome
    let cov_fragment = MidpointCounts::new(mids);
    let rc: Vec<(usize, ReadCount)> = sv_ranges
        .iter()
        .map(|&(id, small, ranges)| {
            let [left_rc, rc, right_rc] = ranges.map(|r| match r {
                None => 0,
                Some((lo, hi)) if small => cov_bases.sum(lo, hi),
                Some((lo, hi)) => cov_fragment.sum(lo, hi),
            });
            (id, ReadCount { left_rc, rc, right_rc })
        })
        .collect();
    Some(ChromResult { events, rc })
}

/// Per-sample genotyping state
struct SampleCounts {
    jct: Vec<JunctionCount>,
    span: Vec<SpanningCount>,
    rc: Vec<ReadCount>,
    dump: Vec<u8>,
}

/// Align one SV in read order. A wave never contains more reads than the remaining
/// capacity, so even parallel waves do no alignment that a reached cap would discard.
fn align_group<F>(mut jct: JunctionCount, mut refc: u32, jobs: Vec<(usize, AlignJob)>, cap: u32,
    parallel_reads: bool, align: &F) -> (JunctionCount, u32, Vec<(usize, String)>)
where F: Fn(AlignJob) -> AlignResult + Sync {
    let mut jobs = jobs.into_iter();
    let mut dumps = Vec::new();
    while jobs.len() > 0 && jct.ref_.len() + jct.alt.len() < cap as usize {
        let remaining = cap as usize - jct.ref_.len() - jct.alt.len();
        let wave: Vec<_> = jobs.by_ref().take(remaining.min(128)).collect();
        let results: Vec<_> = if parallel_reads && wave.len() >= 16 {
            wave.into_par_iter().with_min_len(4).map(|(index, job)| (index, align(job))).collect()
        } else {
            wave.into_iter().map(|(index, job)| (index, align(job))).collect()
        };
        for (index, result) in results {
            if result.kind == b'A' {
                jct.alt.push(result.qual);
                if let Some(dump) = result.dump { dumps.push((index, dump)); }
            } else if result.kind == b'R' {
                refc = refc.wrapping_add(1);
                if refc % 2 == 1 { jct.ref_.push(result.qual); }
            }
        }
    }
    (jct, refc, dumps)
}

struct SampleReplay {
    counts: SampleCounts,
    ref_aligned_read: Vec<u32>,
    ref_aligned_span: Vec<u32>,
    cliptra: HashMap<u64, (bool, u8)>,
}

impl SampleReplay {
    fn new(n: usize) -> Self {
        Self {
            counts: SampleCounts { jct: vec![JunctionCount::default(); n], span: vec![SpanningCount::default(); n],
                rc: vec![ReadCount::default(); n], dump: Vec::new() },
            ref_aligned_read: vec![0; n], ref_aligned_span: vec![0; n], cliptra: HashMap::new(),
        }
    }

    fn full(&self, id: u32, cap: u32) -> bool {
        let j = &self.counts.jct[id as usize];
        j.ref_.len() + j.alt.len() >= cap as usize
    }

    fn process(&mut self, c: &Config, pr: &Probes, mut events: Vec<Event>) {
        let mut groups: BTreeMap<u32, Vec<(usize, AlignJob)>> = BTreeMap::new();
        for (index, ev) in events.iter_mut().enumerate() {
            if matches!(ev, Event::Job(_)) {
                if let Event::Job(job) = std::mem::replace(ev, Event::JunctionDump(None)) {
                    groups.entry(job.id).or_default().push((index, job));
                }
            }
        }
        // Grouping changes execution order only: SV counts are independent, and each
        // group's ordered results return to the original event slots for dump output.
        let parallel_reads = groups.len() < rayon::current_num_threads();
        let work: Vec<_> = groups.into_iter().map(|(id, jobs)| {
            (id, jobs, std::mem::take(&mut self.counts.jct[id as usize]), self.ref_aligned_read[id as usize])
        }).collect();
        let aligned: Vec<_> = work.into_par_iter().map(|(id, jobs, jct, refc)| {
            let (jct, refc, dumps) = align_group(jct, refc, jobs, c.max_geno_read_count,
                parallel_reads, &|job| align_job(c, pr, job));
            (id, jct, refc, dumps)
        }).collect();
        for (id, jct, refc, dumps) in aligned {
            self.counts.jct[id as usize] = jct;
            self.ref_aligned_read[id as usize] = refc;
            for (index, dump) in dumps { events[index] = Event::JunctionDump(Some(dump)); }
        }
        for ev in events {
            match ev {
                Event::Job(_) => unreachable!("all jobs were aligned or capped"),
                Event::JunctionDump(dump) => {
                    if let Some(dump) = dump { self.counts.dump.extend_from_slice(dump.as_bytes()); }
                }
                Event::SpanRef { id, q } => {
                    let id = id as usize;
                    self.ref_aligned_span[id] = self.ref_aligned_span[id].wrapping_add(1);
                    if self.ref_aligned_span[id] % 2 == 1 { self.counts.span[id].ref_.push(q); }
                }
                Event::SpanAlt { id, q, dump } => {
                    if let Some(dump) = dump { self.counts.dump.extend_from_slice(dump.as_bytes()); }
                    self.counts.span[id as usize].alt.push(q);
                }
                Event::TraFirst { hv, clip, q } => { self.cliptra.insert(hv, (clip, q)); }
                Event::TraSecond { hv, q, alts } => {
                    let pq = match self.cliptra.remove(&hv) {
                        None => continue,
                        Some((_, mq)) => mq.min(q),
                    };
                    if (pq as u16) < c.min_geno_qual { continue; }
                    for (id, dump) in alts {
                        if let Some(dump) = dump { self.counts.dump.extend_from_slice(dump.as_bytes()); }
                        self.counts.span[id as usize].alt.push(pq);
                    }
                }
            }
        }
    }
}

/// Sequential admission/replay preserves C++ chromosome and batch boundaries.
fn replay_sample(c: &Config, pr: &Probes, n: usize, chroms: Vec<ChromResult>, batch_size: usize) -> SampleCounts {
    let mut state = SampleReplay::new(n);
    // Dump-mode C++ merges immediately. Per-SV processing reproduces that prefix
    // exactly; bounded chunks only schedule work and preserve all SR/PE interleaving.
    let batch_size = if c.has_dump_file { 4096 } else { batch_size.max(1) };
    for ch in chroms {
        let mut pending = Vec::new();
        let mut jobs = 0;
        for ev in ch.events {
            if let Event::Job(job) = &ev {
                if state.full(job.id, c.max_geno_read_count) { continue; }
                jobs += 1;
            }
            pending.push(ev);
            if jobs >= batch_size {
                state.process(c, pr, std::mem::take(&mut pending));
                jobs = 0;
            }
        }
        state.process(c, pr, pending);
        for (id, r) in ch.rc { state.counts.rc[id] = r; }
    }
    state.counts
}

pub struct Annotation {
    pub rc: Vec<Vec<ReadCount>>,
    pub jct: Vec<Vec<JunctionCount>>,
    pub span: Vec<Vec<SpanningCount>>,
}

/// annotateCoverage
pub fn annotate_coverage(c: &Config, sample_lib: &[LibraryInfo], svs: &mut [SVRecord], hdr: &Header) -> std::io::Result<Annotation> {
    // Test hook: DELLY_BATCH_BASE overrides the 131072 batch base (to check batch-boundary semantics)
    let base: usize = std::env::var("DELLY_BATCH_BASE").ok().and_then(|v| v.parse().ok()).unwrap_or(131072);
    let batch_size = base * c.max_threads as usize;
    let pr = generate_probes(c, hdr, svs);
    log("SV annotation");
    let svs_ro: &[SVRecord] = svs;
    let nf = c.files.len();
    let tasks: Vec<(usize, i32)> = (0..nf).flat_map(|f| (0..hdr.n_targets()).filter(|&r| pr.sv_on_chr[r as usize]).map(move |r| (f, r))).collect();
    let mut chroms: Vec<Option<ChromResult>> = tasks.par_iter().map(|&(f, r)| annotate_chrom(c, &sample_lib[f], svs_ro, &pr, hdr, f, r)).collect();
    let mut per_file: Vec<Vec<ChromResult>> = (0..nf).map(|_| Vec::new()).collect();
    for (k, &(f, _)) in tasks.iter().enumerate() {
        if let Some(r) = chroms[k].take() {
            per_file[f].push(r);
        }
    }
    let per: Vec<SampleCounts> = per_file.into_par_iter().map(|chs| replay_sample(c, &pr, svs_ro.len(), chs, batch_size)).collect();
    if c.has_dump_file {
        let mut all: Vec<u8> = b"#svid\tbam\tqname\tchr\tpos\tmatechr\tmatepos\tmapq\ttype\n".to_vec();
        for s in &per {
            all.extend_from_slice(&s.dump);
        }
        crate::gz::write_gzip_file(&c.dumpfile, &all)?;
    }
    let mut a = Annotation { rc: Vec::new(), jct: Vec::new(), span: Vec::new() };
    for s in per {
        a.rc.push(s.rc);
        a.jct.push(s.jct);
        a.span.push(s.span);
    }
    Ok(a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn job(id: u32, kind: u8, dump: &str) -> AlignJob {
        AlignJob { id, bp_point: 0, svt: 2, sequence: Arc::from(vec![b'A'; 26]),
            qual: kind, dump: Some(Arc::new(DumpRead { file: 0,
                header: Arc::new(Header { names: vec!["chr1".into()], lens: vec![1000], text: String::new() }),
                qname: dump.trim_end().into(), tid: 0, pos: 1, mtid: 0, mpos: 2, qual: 60 })) }
    }

    #[test]
    fn caps_skip_alignment_and_preserve_reference_parity() {
        let kinds: Vec<u8> = b"RRNRARANRRRA".iter().copied().cycle().take(200).collect();
        for cap in [0, 1, 3, 25, 250] {
            for initial_refc in [0u32, 1] {
                let mut expected = JunctionCount::default();
                let mut parity = initial_refc;
                let mut needed = 0;
                let mut dumps = Vec::new();
                for (index, &kind) in kinds.iter().enumerate() {
                    if expected.ref_.len() + expected.alt.len() >= cap as usize { break; }
                    needed += 1;
                    if kind == b'A' {
                        expected.alt.push(25);
                        dumps.push((index, index.to_string()));
                    } else if kind == b'R' {
                        parity += 1;
                        if parity % 2 == 1 { expected.ref_.push(25); }
                    }
                }
                for threads in [1, 4] {
                    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
                    let calls = AtomicUsize::new(0);
                    let jobs = kinds.iter().enumerate().map(|(i, &k)| (i, job(0, k, &i.to_string()))).collect();
                    let (got, refc, got_dumps) = pool.install(|| align_group(JunctionCount::default(), initial_refc,
                        jobs, cap, true, &|j| {
                            calls.fetch_add(1, Ordering::Relaxed);
                            AlignResult { kind: j.qual, qual: 25, dump: j.dump.map(|read| read.qname.clone()) }
                        }));
                    assert_eq!(got.ref_, expected.ref_);
                    assert_eq!(got.alt, expected.alt);
                    assert_eq!(refc, parity);
                    assert_eq!(got_dumps, dumps);
                    assert_eq!(calls.load(Ordering::Relaxed), needed, "cap {}", cap);
                }
            }
        }
    }

    #[test]
    fn replay_retains_dump_order_caps_and_cross_chromosome_mates() {
        let mut pr = Probes { ref_probe: [vec![vec![b'C'; 26]; 2], vec![Vec::new(); 2]],
            cons_probe: [vec![vec![b'A'; 26]; 2], vec![Vec::new(); 2]], bp_region: Vec::new(), sv_on_chr: Vec::new() };
        // Also exercise a rejected alignment, which must not consume capacity.
        pr.ref_probe[1] = vec![vec![b'C'; 26]; 2];
        pr.cons_probe[1] = vec![vec![b'C'; 26]; 2];
        for dump in [false, true] {
            for batch in [1, 2, 100] {
                let mut rejected = job(1, 60, "BAD\n"); rejected.bp_point = 1;
                let chs = vec![ChromResult { events: vec![
                    Event::TraFirst { hv: 7, clip: false, q: 60 },
                    Event::Job(rejected), Event::Job(job(1, 60, "SR1\n")),
                    Event::SpanAlt { id: 0, q: 30, dump: Some("PE\n".into()) },
                    Event::Job(job(0, 60, "SR0\n")), Event::Job(job(1, 60, "CAPPED\n")),
                    Event::SpanRef { id: 0, q: 60 },
                ], rc: vec![(0, ReadCount { rc: 123, ..Default::default() })] },
                ChromResult { events: vec![
                    Event::SpanRef { id: 0, q: 50 },
                    Event::TraSecond { hv: 7, q: 40, alts: vec![(1, Some("TRA\n".into()))] },
                    Event::Job(job(0, 60, "CAPPED2\n")),
                ], rc: Vec::new() }];
                let c = Config { has_dump_file: dump, max_geno_read_count: 1, files: vec!["test.bam".into()], ..Default::default() };
                let got = replay_sample(&c, &pr, 2, chs, batch);
                assert_eq!(got.jct[0].alt, vec![25]);
                assert_eq!(got.jct[1].alt, vec![25]);
                assert_eq!(got.span[0].ref_, vec![60]);
                assert_eq!(got.span[0].alt, vec![30]);
                assert_eq!(got.span[1].alt, vec![40]);
                assert_eq!(got.rc[0].rc, 123);
                let text = String::from_utf8(got.dump).unwrap();
                let labels: Vec<_> = text.lines().map(|line| line.split('\t').nth(2).unwrap_or(line)).collect();
                assert_eq!(labels, ["SR1", "PE", "SR0", "TRA"]);
            }
        }
    }

    #[test]
    fn compact_coverage_matches_dense_arrays() {
        let len = 5000u32;
        let mut seed = 7u64;
        let mut next = |m: u32| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((seed >> 33) % m as u64) as u32
        };
        // Overlapping, adjacent, empty and clipped ranges
        let mut ranges: Vec<(i32, i32)> = (0..60).map(|_| (next(len + 200) as i32 - 100, next(len + 200) as i32 - 100)).collect();
        ranges.extend([(0, 10), (10, 20), (2400, 2600), (2501, 2502), (4990, 6000), (-5, 3), (30, 30)]);
        let wanted: Vec<(u32, u32)> = ranges.iter().filter_map(|&(a, b)| rc_range(a, b, len)).collect();
        let mut dense = vec![0u16; len as usize];
        let mut windowed = WindowedCounts::new(wanted.clone());
        let mut dense_mid = vec![0u16; len as usize];
        let mut mids = Vec::new();
        // A hot spot exceeds the uint16_t saturation limit
        for i in 0..70_000u32 {
            let (s, l) = if i < 66_000 { (2500, 3) } else { (next(len), 1 + next(150)) };
            let e = (s + l).min(len);
            windowed.add(s, e);
            for v in &mut dense[s as usize..e as usize] {
                if *v < COV_SATURATION { *v += 1; }
            }
            let m = if i < 66_000 { 17 } else { next(len) };
            mids.push(m);
            if dense_mid[m as usize] < COV_SATURATION { dense_mid[m as usize] += 1; }
        }
        let mids = MidpointCounts::new(mids);
        for &(lo, hi) in &wanted {
            let sum = |d: &[u16]| d[lo as usize..hi as usize].iter().fold(0i32, |s, &v| s.wrapping_add(v as i32));
            assert_eq!(windowed.sum(lo, hi), sum(&dense), "bases {}..{}", lo, hi);
            assert_eq!(mids.sum(lo, hi), sum(&dense_mid), "midpoints {}..{}", lo, hi);
        }
    }

    #[test]
    fn reference_reads_do_not_format_dump_mate_names() {
        let pr = Probes { ref_probe: [vec![vec![b'A'; 26]], vec![Vec::new()]],
            cons_probe: [vec![vec![b'C'; 26]], vec![Vec::new()]], bp_region: Vec::new(), sv_on_chr: Vec::new() };
        let mut candidate = job(0, 60, "unpaired");
        Arc::get_mut(candidate.dump.as_mut().unwrap()).unwrap().mtid = -1;
        let result = align_job(&Config::default(), &pr, candidate);
        assert_eq!(result.kind, b'R');
        assert!(result.dump.is_none());
    }
}
