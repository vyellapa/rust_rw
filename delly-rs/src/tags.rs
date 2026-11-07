//! Core records and helpers (port of tags.h and parts of util.h / cluster.h).

use crate::cxx_sort::CppLess;
use crate::hts::*;

pub const DELLY_SVT_TRANS: i32 = 5;

#[inline]
pub fn is_translocation(svt: i32) -> bool {
    DELLY_SVT_TRANS <= svt && svt < 9
}

#[inline]
pub fn rec_translocation(rec: &Record) -> bool {
    rec.tid() != rec.mtid()
}

#[inline]
pub fn span_orientation(svt: i32) -> u8 {
    if is_translocation(svt) {
        (svt - DELLY_SVT_TRANS) as u8
    } else {
        svt as u8
    }
}

/// Junction record
#[derive(Clone, Copy, Debug)]
pub struct Junction {
    pub forward: bool,
    pub scleft: bool,
    pub refidx: i32,
    pub rstart: i32,
    pub refpos: i32,
    pub seqpos: i32,
    pub qual: u16,
}

impl CppLess for Junction {
    #[inline]
    fn lt(&self, j2: &Self) -> bool {
        if self.seqpos != j2.seqpos {
            return self.seqpos < j2.seqpos;
        }
        if self.refidx != j2.refidx {
            return self.refidx < j2.refidx;
        }
        if self.refpos != j2.refpos {
            return self.refpos < j2.refpos;
        }
        if self.scleft != j2.scleft {
            return !self.scleft & j2.scleft;
        }
        if self.rstart != j2.rstart {
            return self.rstart < j2.rstart;
        }
        if self.forward != j2.forward {
            return !self.forward & j2.forward;
        }
        self.qual < j2.qual
    }
}

/// Split-read record
#[derive(Clone, Copy, Debug)]
pub struct SRBamRecord {
    pub chr: i32,
    pub pos: i32,
    pub chr2: i32,
    pub pos2: i32,
    pub rstart: i32,
    pub sstart: i32,
    pub qual: i32,
    pub inslen: i32,
    pub svid: i32,
    pub primary_chr: i32,
    pub id: u64,
}

impl SRBamRecord {
    #[allow(clippy::too_many_arguments)]
    pub fn new(c: i32, p: i32, c2: i32, p2: i32, rst: i32, sst: i32, qval: i32, il: i32, idval: u64) -> Self {
        SRBamRecord { chr: c, pos: p, chr2: c2, pos2: p2, rstart: rst, sstart: sst, qual: qval, inslen: il, svid: -1, primary_chr: -1, id: idval }
    }
}

impl CppLess for SRBamRecord {
    #[inline]
    fn lt(&self, s: &Self) -> bool {
        if self.chr != s.chr {
            return self.chr < s.chr;
        }
        if self.pos != s.pos {
            return self.pos < s.pos;
        }
        if self.chr2 != s.chr2 {
            return self.chr2 < s.chr2;
        }
        if self.pos2 != s.pos2 {
            return self.pos2 < s.pos2;
        }
        if self.rstart != s.rstart {
            return self.rstart < s.rstart;
        }
        if self.sstart != s.sstart {
            return self.sstart < s.sstart;
        }
        if self.inslen != s.inslen {
            return self.inslen < s.inslen;
        }
        if self.qual != s.qual {
            return self.qual < s.qual;
        }
        self.id < s.id
    }
}

#[derive(Clone, Debug, Default)]
pub struct SVAnno {
    pub is_rc: bool,
    pub seq_type: i32,
    pub hom_len: i32,
    pub tr_period: i32,
    pub tr_copies: f32,
    pub from_sites: bool,
}

/// Structural variant record
#[derive(Clone, Debug)]
pub struct SVRecord {
    pub chr: i32,
    pub sv_start: i32,
    pub chr2: i32,
    pub sv_end: i32,
    pub ciposlow: i32,
    pub ciposhigh: i32,
    pub ciendlow: i32,
    pub ciendhigh: i32,
    pub sr_support: i32,
    pub sr_map_quality: i32,
    pub mapq: i32,
    pub ins_len: i32,
    pub svt: i32,
    pub id: i32,
    pub hom_len: i32,
    pub pe_support: i32,
    pub pe_map_quality: i32,
    pub cons_bp: i32,
    pub alleleid: i32,
    pub nallele: i32,
    pub tr_start: i32,
    pub tr_end: i32,
    pub sr_align_quality: f32,
    pub precise: bool,
    pub alleles: String,
    pub consensus: String,
    pub anno: SVAnno,
}

impl Default for SVRecord {
    fn default() -> Self {
        SVRecord {
            chr: 0,
            sv_start: 0,
            chr2: 0,
            sv_end: 0,
            ciposlow: 0,
            ciposhigh: 0,
            ciendlow: 0,
            ciendhigh: 0,
            sr_support: 0,
            sr_map_quality: 0,
            mapq: 0,
            ins_len: 0,
            svt: -1,
            id: 0,
            hom_len: 0,
            pe_support: 0,
            pe_map_quality: 0,
            cons_bp: 0,
            alleleid: -1,
            nallele: 1,
            tr_start: -1,
            tr_end: -1,
            sr_align_quality: 0.0,
            precise: false,
            alleles: String::new(),
            consensus: String::new(),
            anno: SVAnno::default(),
        }
    }
}

impl SVRecord {
    /// StructuralVariantRecord(c, s, e): lookup key
    pub fn key(c: i32, s: i32, e: i32) -> Self {
        SVRecord { chr: c, sv_start: s, chr2: c, sv_end: e, ..Default::default() }
    }

    /// The split-read constructor (precise = true)
    #[allow(clippy::too_many_arguments)]
    pub fn new_sr(c1: i32, s: i32, c2: i32, e: i32, cipl: i32, ciph: i32, ciel: i32, cieh: i32, sup: i32, srmapq: i32, qval: i32, ilen: i32, svtype: i32, idval: i32) -> Self {
        SVRecord {
            chr: c1,
            sv_start: s,
            chr2: c2,
            sv_end: e,
            ciposlow: cipl,
            ciposhigh: ciph,
            ciendlow: ciel,
            ciendhigh: cieh,
            sr_support: sup,
            sr_map_quality: srmapq,
            mapq: qval,
            ins_len: ilen,
            svt: svtype,
            id: idval,
            precise: true,
            ..Default::default()
        }
    }
}

impl CppLess for SVRecord {
    #[inline]
    fn lt(&self, s: &Self) -> bool {
        if self.chr != s.chr {
            return self.chr < s.chr;
        }
        if self.sv_start != s.sv_start {
            return self.sv_start < s.sv_start;
        }
        if self.chr2 != s.chr2 {
            return self.chr2 < s.chr2;
        }
        if self.sv_end != s.sv_end {
            return self.sv_end < s.sv_end;
        }
        if self.pe_support != s.pe_support {
            return self.pe_support > s.pe_support;
        }
        if self.sr_support != s.sr_support {
            return self.sr_support > s.sr_support;
        }
        if self.svt != s.svt {
            return self.svt < s.svt;
        }
        if self.ins_len != s.ins_len {
            return self.ins_len < s.ins_len;
        }
        if self.mapq != s.mapq {
            return self.mapq > s.mapq;
        }
        if self.sr_map_quality != s.sr_map_quality {
            return self.sr_map_quality > s.sr_map_quality;
        }
        // std::string operator< (lexicographic, unsigned bytes via char_traits<char>::compare = memcmp)
        self.consensus.as_bytes() < s.consensus.as_bytes()
    }
}

/// Breakpoint
#[allow(dead_code)]
#[derive(Clone, Debug, Default)]
pub struct Breakpoint {
    pub sv_start_beg: i32,
    pub sv_start_end: i32,
    pub sv_end_beg: i32,
    pub sv_end_end: i32,
    pub sv_start: i32,
    pub sv_end: i32,
    pub pe_support: i32,
    pub svt: i32,
    pub chr: i32,
    pub chr2: i32,
    pub part1: Vec<u8>,
}

impl Breakpoint {
    pub fn from_sv(sv: &SVRecord) -> Self {
        Breakpoint {
            sv_start_beg: sv.sv_start,
            sv_start_end: sv.sv_start,
            sv_end_beg: sv.sv_end,
            sv_end_end: sv.sv_end,
            sv_start: sv.sv_start,
            sv_end: sv.sv_end,
            pe_support: sv.pe_support,
            svt: sv.svt,
            chr: sv.chr,
            chr2: sv.chr2,
            part1: Vec::new(),
        }
    }
}

/// _initBreakpoint
pub fn init_breakpoint(target_len: &[u32], bp: &mut Breakpoint, boundary: i32, svt: i32) {
    if is_translocation(svt) || svt == 4 {
        bp.sv_start_beg = 0.max(bp.sv_start - boundary);
        bp.sv_start_end = (target_len[bp.chr as usize] as i32).min(bp.sv_start + boundary);
        bp.sv_end_beg = 0.max(bp.sv_end - boundary);
        bp.sv_end_end = (target_len[bp.chr2 as usize] as i32).min(bp.sv_end + boundary);
    } else {
        bp.sv_start_beg = 0.max(bp.sv_start - boundary);
        bp.sv_start_end = (bp.sv_start + boundary).min((bp.sv_start + bp.sv_end) / 2);
        bp.sv_end_beg = ((bp.sv_start + bp.sv_end) / 2 + 1).max(bp.sv_end - boundary);
        bp.sv_end_end = (target_len[bp.chr2 as usize] as i32).min(bp.sv_end + boundary);
    }
}

#[inline]
pub fn min_coord(position: i32, mate_position: i32, svt: i32) -> i32 {
    if is_translocation(svt) {
        position
    } else {
        position.min(mate_position)
    }
}

#[inline]
pub fn max_coord(position: i32, mate_position: i32, svt: i32) -> i32 {
    if is_translocation(svt) {
        mate_position
    } else {
        position.max(mate_position)
    }
}

/// Short-read SV size check
#[inline]
pub fn sv_size_check(s: i32, e: i32, svt: i32) -> bool {
    match svt {
        0 | 1 | 2 => (e - s) >= 300,
        3 => (e - s) >= 100,
        _ => true,
    }
}

/// Long-read SV size check (used by the split-read clique search)
#[inline]
pub fn sv_size_check_ins(s: i32, e: i32, svt: i32, inslen: i32) -> bool {
    match svt {
        0 | 1 | 2 | 3 => (e - s) >= 15,
        4 => inslen >= 15,
        _ => true,
    }
}

/// getSVType
#[inline]
pub fn get_sv_type(rec: &Record) -> u8 {
    let f = rec.flag();
    if f & BAM_FREVERSE == 0 {
        if f & BAM_FMREVERSE == 0 {
            0
        } else if rec.pos() < rec.mpos() {
            2
        } else {
            3
        }
    } else if f & BAM_FMREVERSE == 0 {
        if rec.pos() > rec.mpos() {
            2
        } else {
            3
        }
    } else {
        1
    }
}

/// _isizeMappingPos
pub fn isize_mapping_pos(rec: &Record, isize: i32) -> i32 {
    if rec_translocation(rec) {
        let orient = get_sv_type(rec);
        if orient == 0 {
            DELLY_SVT_TRANS
        } else if orient == 1 {
            DELLY_SVT_TRANS + 1
        } else if rec.tid() > rec.mtid() {
            if rec.flag() & BAM_FREVERSE == 0 {
                DELLY_SVT_TRANS + 2
            } else {
                DELLY_SVT_TRANS + 3
            }
        } else if rec.flag() & BAM_FREVERSE == 0 {
            DELLY_SVT_TRANS + 3
        } else {
            DELLY_SVT_TRANS + 2
        }
    } else {
        if rec.pos() == rec.mpos() {
            return -1;
        }
        let orient = get_sv_type(rec);
        if orient == 0 {
            0
        } else if orient == 1 {
            1
        } else if orient == 2 {
            // isize > std::abs(rec->core.isize): int32 vs hts_pos_t (int64) comparison
            if (isize as i64) > rec.isize().abs() {
                -1
            } else {
                2
            }
        } else {
            if (rec.pos() - rec.mpos()).abs() < 100 {
                return -1;
            }
            3
        }
    }
}

/// _firstPairObs
#[inline]
pub fn first_pair_obs(rec: &Record, last_aligned_pos_reads: &std::collections::HashSet<u64>) -> bool {
    if rec.tid() == rec.mtid() {
        (rec.pos() < rec.mpos()) || ((rec.pos() == rec.mpos()) && !last_aligned_pos_reads.contains(&(crate::cxx_hash::hash_string(rec.qname()) as u64)))
    } else {
        rec.tid() < rec.mtid()
    }
}

/// _pairsDisagree
#[allow(clippy::too_many_arguments)]
pub fn pairs_disagree(p1min: i32, p1max: i32, p1rl: i32, p1mis: i32, p2min: i32, p2max: i32, p2rl: i32, p2mis: i32, svt: i32) -> bool {
    if is_translocation(svt) {
        let ct = span_orientation(svt);
        if ct % 2 == 0 {
            if (p2min + p2rl - p1min) > p1mis {
                return true;
            }
            if ct >= 2 {
                if p2max < p1max {
                    if (p1max + p1rl - p2max) > p1mis {
                        return true;
                    }
                } else if (p2max + p2rl - p1max) > p2mis {
                    return true;
                }
            } else if p2max < p1max {
                if (p1max + p1rl - p2max) > p2mis {
                    return true;
                }
            } else if (p2max + p2rl - p1max) > p1mis {
                return true;
            }
        } else {
            if (p2min + p2rl - p1min) > p2mis {
                return true;
            }
            if ct >= 2 {
                if p2max < p1max {
                    if (p1max + p1rl - p2max) > p2mis {
                        return true;
                    }
                } else if (p2max + p2rl - p1max) > p1mis {
                    return true;
                }
            } else if p2max < p1max {
                if (p1max + p1rl - p2max) > p1mis {
                    return true;
                }
            } else if (p2max + p2rl - p1max) > p2mis {
                return true;
            }
        }
        false
    } else if svt < 2 {
        if svt == 0 {
            if (p2min + p2rl - p1min) > p1mis {
                return true;
            }
            if (p2max < p1max) && ((p1max + p1rl - p2max) > p2mis) {
                return true;
            }
            if (p2max >= p1max) && ((p2max + p2rl - p1max) > p1mis) {
                return true;
            }
        } else {
            if (p2min + p2rl - p1min) > p2mis {
                return true;
            }
            if (p2max < p1max) && ((p1max + p1rl - p2max) > p1mis) {
                return true;
            }
            if (p2max >= p1max) && ((p2max + p2rl - p1max) > p2mis) {
                return true;
            }
        }
        false
    } else if svt == 2 {
        if (p2min + p2rl - p1min) > p1mis {
            return true;
        }
        if (p2max < p1max) && ((p1max + p1rl - p2max) > p1mis) {
            return true;
        }
        if (p2max >= p1max) && ((p2max + p2rl - p1max) > p2mis) {
            return true;
        }
        if (p1max < p2min) || (p2max < p1min) {
            return true;
        }
        false
    } else if svt == 3 {
        if (p2min + p2rl - p1min) > p2mis {
            return true;
        }
        if (p2max < p1max) && ((p1max + p1rl - p2max) > p2mis) {
            return true;
        }
        if (p2max >= p1max) && ((p2max + p2rl - p1max) > p1mis) {
            return true;
        }
        false
    } else {
        false
    }
}

/// Library information
#[derive(Clone, Debug, Default)]
pub struct LibraryInfo {
    pub rs: i32,
    pub median: i32,
    pub mad: i32,
    pub min_normal_isize: i32,
    pub min_isize_cutoff: i32,
    pub max_normal_isize: i32,
    pub max_isize_cutoff: i32,
    pub abnormal_pairs: u32,
}

/// Reduced paired-end alignment record (cluster.h BamAlignRecord)
#[allow(dead_code)]
#[derive(Clone, Copy, Debug)]
pub struct BamAlignRecord {
    pub tid: i32,
    pub pos: i32,
    pub mtid: i32,
    pub mpos: i32,
    pub alen: i32,
    pub malen: i32,
    pub median: i32,
    pub mad: i32,
    pub max_normal_isize: i32,
    pub flag: u32,
    pub map_quality: u8,
}

impl CppLess for BamAlignRecord {
    #[inline]
    fn lt(&self, s2: &Self) -> bool {
        if self.tid == self.mtid {
            let lo1 = self.pos.min(self.mpos);
            let lo2 = s2.pos.min(s2.mpos);
            let hi1 = self.pos.max(self.mpos);
            let hi2 = s2.pos.max(s2.mpos);
            if lo1 != lo2 {
                return lo1 < lo2;
            }
            if hi1 != hi2 {
                return hi1 < hi2;
            }
            if self.max_normal_isize != s2.max_normal_isize {
                return self.max_normal_isize < s2.max_normal_isize;
            }
        } else {
            if self.pos != s2.pos {
                return self.pos < s2.pos;
            }
            if self.mpos != s2.mpos {
                return self.mpos < s2.mpos;
            }
            if self.max_normal_isize != s2.max_normal_isize {
                return self.max_normal_isize < s2.max_normal_isize;
            }
        }
        if self.tid != s2.tid {
            return self.tid < s2.tid;
        }
        if self.mtid != s2.mtid {
            return self.mtid < s2.mtid;
        }
        if self.pos != s2.pos {
            return self.pos < s2.pos;
        }
        if self.flag != s2.flag {
            return self.flag < s2.flag;
        }
        if self.map_quality != s2.map_quality {
            return self.map_quality > s2.map_quality;
        }
        if self.alen != s2.alen {
            return self.alen < s2.alen;
        }
        self.malen < s2.malen
    }
}

/// _addID
pub fn add_id(svt: i32) -> &'static str {
    match svt {
        0 | 1 => "INV",
        2 => "DEL",
        3 => "DUP",
        4 => "INS",
        9 => "CNV",
        _ => "BND",
    }
}

/// _addOrientation
pub fn add_orientation(svt: i32) -> &'static str {
    match span_orientation(svt) {
        0 => "3to3",
        1 => "5to5",
        2 => "3to5",
        3 => "5to3",
        _ => "NtoN",
    }
}

/// _decodeOrientation(value, svt)
pub fn decode_orientation(value: &str, svt: &str) -> i32 {
    if svt == "BND" {
        match value {
            "3to3" => DELLY_SVT_TRANS,
            "5to5" => DELLY_SVT_TRANS + 1,
            "3to5" => DELLY_SVT_TRANS + 2,
            "5to3" => DELLY_SVT_TRANS + 3,
            _ => -1,
        }
    } else if svt == "CNV" {
        9
    } else {
        match value {
            "3to3" => 0,
            "5to5" => 1,
            "3to5" => 2,
            "5to3" => 3,
            _ => 4,
        }
    }
}
