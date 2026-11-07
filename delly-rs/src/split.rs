//! Port of split.h: consensus-to-reference split alignment and breakpoint refinement.

use crate::align::{long_needle, longest_homology, Align, AlignConfig};
use crate::config::{Config, DnaScore};
use crate::edlib::{self, Mode, Task, EDLIB_EDOP_DELETE, EDLIB_EDOP_INSERT};
use crate::hts::RefSeq;
use crate::tags::*;
use crate::util::{reverse_complement, to_upper};

#[derive(Clone, Copy, Debug, Default)]
pub struct AlignDescriptor {
    pub c_start: i32,
    pub c_end: i32,
    pub r_start: i32,
    pub r_end: i32,
    pub hom_left: i32,
    pub hom_right: i32,
    pub perc_id: f32,
}

/// std::string::substr(pos) / substr(pos, n) semantics (pos <= len assumed)
#[inline]
pub fn substr(s: &[u8], pos: usize, n: usize) -> &[u8] {
    let end = if n > s.len() - pos { s.len() } else { pos + n };
    &s[pos..end]
}

#[inline]
pub fn substr_from(s: &[u8], pos: usize) -> &[u8] {
    &s[pos..]
}

/// _adjustOrientation
pub fn adjust_orientation(seq: &mut [u8], bp_point: bool, svt: i32) {
    if is_translocation(svt) {
        let ct = span_orientation(svt);
        if (ct == 0 && bp_point) || (ct == 1 && !bp_point) {
            reverse_complement(seq);
        }
    } else if svt == 0 {
        if bp_point {
            reverse_complement(seq);
        }
    } else if svt == 1 && !bp_point {
        reverse_complement(seq);
    }
}

/// The reversed-complement loop used inside _getSVRef (keeps unknown characters of the original position)
fn revcomp_keep(s: &[u8]) -> Vec<u8> {
    let mut out = s.to_vec();
    for (i, &c) in s.iter().rev().enumerate() {
        match c {
            b'A' => out[i] = b'T',
            b'C' => out[i] = b'G',
            b'G' => out[i] = b'C',
            b'T' => out[i] = b'A',
            b'N' => out[i] = b'N',
            _ => {}
        }
    }
    out
}

#[inline]
fn up(r: &RefSeq, a: i32, b: i32) -> Vec<u8> {
    to_upper(&r.seq()[a as usize..b as usize])
}

/// _getSVRef
pub fn get_sv_ref(c: &Config, rf: &RefSeq, bp: &Breakpoint, ref_index: i32, svt: i32) -> Vec<u8> {
    if is_translocation(svt) {
        let ct = span_orientation(svt);
        if bp.chr == ref_index {
            if ct == 0 || ct == 2 {
                let mut v = up(rf, bp.sv_start_beg, bp.sv_start_end);
                v.extend_from_slice(&bp.part1);
                v
            } else if ct == 1 {
                let mut v = revcomp_keep(&up(rf, bp.sv_start_beg, bp.sv_start_end));
                v.extend_from_slice(&bp.part1);
                v
            } else {
                let mut v = bp.part1.clone();
                v.extend_from_slice(&up(rf, bp.sv_start_beg, bp.sv_start_end));
                v
            }
        } else if ct == 0 {
            revcomp_keep(&up(rf, bp.sv_end_beg, bp.sv_end_end))
        } else {
            up(rf, bp.sv_end_beg, bp.sv_end_end)
        }
    } else if svt == 2 {
        if bp.sv_end - bp.sv_start <= c.indelsize {
            up(rf, bp.sv_start_beg, bp.sv_end_end)
        } else {
            let mut v = up(rf, bp.sv_start_beg, bp.sv_start_end);
            v.extend_from_slice(&up(rf, bp.sv_end_beg, bp.sv_end_end));
            v
        }
    } else if svt == 4 {
        up(rf, bp.sv_start_beg, bp.sv_end_end)
    } else if svt == 3 {
        let mut v = up(rf, bp.sv_end_beg, bp.sv_end_end);
        v.extend_from_slice(&up(rf, bp.sv_start_beg, bp.sv_start_end));
        v
    } else if svt == 0 {
        let wide = (bp.sv_end - bp.sv_start) > c.min_cons_window;
        let str_end = if wide { up(rf, bp.sv_end_beg, bp.sv_end_end) } else { up(rf, bp.sv_start, bp.sv_end_end) };
        let rc = revcomp_keep(&str_end);
        let mut v = up(rf, bp.sv_start_beg, bp.sv_start_end);
        v.extend_from_slice(&rc);
        if !wide {
            v.extend_from_slice(&up(rf, bp.sv_end, bp.sv_end_end));
        }
        v
    } else if svt == 1 {
        let wide = (bp.sv_end - bp.sv_start) > c.min_cons_window;
        let str_start = if wide { up(rf, bp.sv_start_beg, bp.sv_start_end) } else { up(rf, bp.sv_start_beg, bp.sv_end) };
        let rc = revcomp_keep(&str_start);
        if wide {
            let mut v = rc;
            v.extend_from_slice(&up(rf, bp.sv_end_beg, bp.sv_end_end));
            v
        } else {
            let mut v = up(rf, bp.sv_start_beg, bp.sv_start);
            v.extend_from_slice(&rc);
            v.extend_from_slice(&up(rf, bp.sv_end_beg, bp.sv_end_end));
            v
        }
    } else {
        Vec::new()
    }
}

/// _coordTransform; results are `unsigned int` in C++ (wrapping)
fn coord_transform(c: &Config, ref_len: usize, sv: &Breakpoint, ad: &AlignDescriptor, svt: i32) -> Option<(u32, u32)> {
    let rl = ref_len as i64;
    let (rs, re) = (ad.r_start as i64, ad.r_end as i64);
    let w = |v: i64| v as u32;
    if is_translocation(svt) {
        let ct = span_orientation(svt);
        match ct {
            0 => {
                let annealed = (sv.sv_start_end - sv.sv_start_beg) as i64;
                if rs >= annealed || re < annealed {
                    return None;
                }
                Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_end_beg as i64 + (rl - re) + 1)))
            }
            1 => {
                let annealed = (sv.sv_start_end - sv.sv_start_beg) as i64;
                if rs >= annealed || re < annealed {
                    return None;
                }
                Some((w(sv.sv_start_beg as i64 + (annealed - rs) + 1), w(sv.sv_end_beg as i64 + (re - annealed))))
            }
            2 => {
                let annealed = (sv.sv_start_end - sv.sv_start_beg) as i64;
                if rs >= annealed || re < annealed {
                    return None;
                }
                Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_end_beg as i64 + (re - annealed))))
            }
            3 => {
                let annealed = (sv.sv_end_end - sv.sv_end_beg) as i64;
                if rs >= annealed || re < annealed {
                    return None;
                }
                Some((w(sv.sv_start_beg as i64 + (re - annealed)), w(sv.sv_end_beg as i64 + rs)))
            }
            _ => None,
        }
    } else if svt == 2 {
        if sv.sv_end - sv.sv_start > c.indelsize {
            let annealed = (sv.sv_start_end - sv.sv_start_beg) as i64;
            if rs >= annealed || re < annealed {
                return None;
            }
            Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_end_beg as i64 + (re - annealed))))
        } else {
            Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_start_beg as i64 + re)))
        }
    } else if svt == 3 {
        let annealed = (sv.sv_end_end - sv.sv_end_beg) as i64;
        if rs >= annealed || re < annealed {
            return None;
        }
        Some((w(sv.sv_start_beg as i64 + (re - annealed)), w(sv.sv_end_beg as i64 + rs)))
    } else if svt == 0 {
        let annealed = (sv.sv_start_end - sv.sv_start_beg) as i64;
        if rs >= annealed || re < annealed {
            return None;
        }
        if (sv.sv_end - sv.sv_start) > c.min_cons_window {
            Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_end_beg as i64 + (rl - re) + 1)))
        } else {
            Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_end_end as i64 - (re - annealed))))
        }
    } else if svt == 1 {
        let annealed = if (sv.sv_end - sv.sv_start) > c.min_cons_window {
            (sv.sv_start_end - sv.sv_start_beg) as i64
        } else {
            ((sv.sv_start - sv.sv_start_beg) + (sv.sv_end - sv.sv_start_beg)) as i64
        };
        if rs >= annealed || re < annealed {
            return None;
        }
        Some((w(sv.sv_start_beg as i64 + (annealed - rs) + 1), w(sv.sv_end_beg as i64 + (re - annealed))))
    } else if svt == 4 {
        Some((w(sv.sv_start_beg as i64 + rs), w(sv.sv_start_beg as i64 + re)))
    } else {
        // returns true without touching the (zero-initialised) outputs
        Some((0, 0))
    }
}

#[inline]
fn valid_sr_alignment(c_start: i32, c_end: i32, r_start: i32, r_end: i32, svt: i32) -> bool {
    if svt == 4 {
        (r_end - r_start) < 5 && (c_end - c_start) > 15
    } else {
        (c_end - c_start) < 5 && (r_end - r_start) > 15
    }
}

#[inline]
fn check_sv_gap(ref_gap: i32, old_ref_gap: i32, var_gap: i32, old_var_gap: i32, svt: i32) -> bool {
    if svt == 4 {
        var_gap > old_var_gap
    } else {
        ref_gap > old_ref_gap
    }
}

fn find_homology(consensus: &[u8], sv_ref: &[u8], ad: &mut AlignDescriptor, svt: i32) {
    if svt == 4 {
        ad.hom_right = longest_homology(substr_from(consensus, ad.c_start as usize), substr_from(sv_ref, (ad.r_end - 1) as usize), -1);
        let mut pre_c = substr(consensus, 0, (ad.c_end - 1) as usize).to_vec();
        let mut pre_r = substr(sv_ref, 0, ad.r_start as usize).to_vec();
        pre_c.reverse();
        pre_r.reverse();
        ad.hom_left = longest_homology(&pre_c, &pre_r, -1);
    } else {
        ad.hom_right = longest_homology(substr_from(consensus, (ad.c_end - 1) as usize), substr_from(sv_ref, ad.r_start as usize), -1);
        let mut pre_c = substr(consensus, 0, ad.c_start as usize).to_vec();
        let mut pre_r = substr(sv_ref, 0, (ad.r_end - 1) as usize).to_vec();
        pre_c.reverse();
        pre_r.reverse();
        ad.hom_left = longest_homology(&pre_c, &pre_r, -1);
    }
}

fn percent_identity(al: &Align, gs: i32, ge: i32) -> f32 {
    let mut var_seen = false;
    let mut ref_seen = false;
    let mut gap_mm: u32 = 0;
    let mut mm: u32 = 0;
    let mut ma: u32 = 0;
    let mut in_gap = false;
    for j in 0..al.cols {
        let ji = j as i32;
        if ji < gs || ji > ge {
            let a0 = al.get(0, j);
            let a1 = al.get(1, j);
            if a0 != b'-' {
                var_seen = true;
            }
            if a1 != b'-' {
                ref_seen = true;
            }
            if a0 == b'-' || a1 == b'-' {
                if ref_seen && var_seen {
                    if !in_gap {
                        in_gap = true;
                        gap_mm = 0;
                    }
                    gap_mm += 1;
                }
            } else {
                if in_gap {
                    mm += gap_mm;
                    in_gap = false;
                }
                if a0 == a1 {
                    ma += 1;
                } else {
                    mm += 1;
                }
            }
        }
    }
    ma as f32 / (ma + mm) as f32
}

/// _findSplit
pub fn find_split(c: &Config, consensus: &[u8], sv_ref: &[u8], al: &Align, ad: &mut AlignDescriptor, svt: i32) -> bool {
    let mut gs: i32 = 0;
    let mut ge: i32 = 0;
    let mut ref_index: i32 = 0;
    let mut var_index: i32 = 0;
    let mut gap_start_ref: i32 = 0;
    let mut gap_start_var: i32 = 0;
    let mut a1: i32 = 0;
    let mut in_gap = false;
    for j in 0..al.cols {
        let x0 = al.get(0, j);
        let x1 = al.get(1, j);
        if x0 != b'-' {
            var_index += 1;
        }
        if x1 != b'-' {
            ref_index += 1;
        }
        if (x0 == b'-' || x1 == b'-') && ref_index > 0 && var_index > 0 {
            if !in_gap {
                gap_start_var = if x0 != b'-' { var_index - 1 } else { var_index };
                gap_start_ref = if x1 != b'-' { ref_index - 1 } else { ref_index };
                a1 = j as i32;
                in_gap = true;
            }
        } else {
            if in_gap && check_sv_gap(ref_index - gap_start_ref, ad.r_end - ad.r_start, var_index - gap_start_var, ad.c_end - ad.c_start, svt) {
                ad.r_start = gap_start_ref;
                ad.r_end = ref_index;
                ad.c_start = gap_start_var;
                ad.c_end = var_index;
                gs = a1;
                ge = j as i32 - 1;
            }
            in_gap = false;
        }
    }
    if ad.r_end <= ad.r_start {
        return false;
    }
    if !valid_sr_alignment(ad.c_start, ad.c_end, ad.r_start, ad.r_end, svt) {
        return false;
    }
    ad.perc_id = percent_identity(al, gs, ge);
    if ad.perc_id < c.flank_quality {
        return false;
    }
    find_homology(consensus, sv_ref, ad, svt);
    if (ad.hom_left + c.minimum_flank_size > ad.c_start) || (var_index < ad.c_end + ad.hom_right + c.minimum_flank_size) {
        return false;
    }
    if (ad.hom_left + c.minimum_flank_size > ad.r_start) || (ref_index < ad.r_end + ad.hom_right + c.minimum_flank_size) {
        return false;
    }
    true
}

/// editDistanceVec
fn edit_distance_vec(seq_i: &[u8], seq_j: &[u8], cig: &edlib::AlignResult) -> Vec<u32> {
    let mut dist = vec![0u32; seq_i.len()];
    let mut t_idx: i32 = -1;
    let mut q_idx: i32 = -1;
    let mut ed: u32 = 0;
    for &op in &cig.alignment {
        if op == EDLIB_EDOP_DELETE {
            t_idx += 1;
            ed += 1;
        } else if op == EDLIB_EDOP_INSERT {
            q_idx += 1;
            ed += 1;
            dist[q_idx as usize] = ed;
        } else {
            t_idx += 1;
            q_idx += 1;
            if seq_i[q_idx as usize] != seq_j[t_idx as usize] {
                ed += 1;
            }
            dist[q_idx as usize] = ed;
        }
    }
    dist
}

/// glueAlignment (HW mode)
fn glue_alignment(query: &[u8], target: &[u8], gaplen: u32, left: &edlib::AlignResult, right: &edlib::AlignResult) -> Align {
    let mut t_idx: i32 = left.end_location;
    let mut q_idx: i32 = -1;
    let mut missing_start: u32 = 0;
    for &op in &left.alignment {
        if op != EDLIB_EDOP_INSERT {
            t_idx -= 1;
        }
    }
    if t_idx >= 0 {
        missing_start = (t_idx + 1) as u32;
    }
    let mut missing_end: u32 = right.end_location as u32;
    if (missing_end as usize) < target.len() {
        missing_end = (target.len() as u32).wrapping_sub(missing_end).wrapping_sub(1);
    }
    let ll = left.alignment.len();
    let rl = right.alignment.len();
    let ms = missing_start as usize;
    let gl = gaplen as usize;
    let me = missing_end as usize;
    let mut al = Align::new(2, ms + ll + gl + rl + me);
    for j in 0..ms {
        al.set(1, j, target[j]);
        al.set(0, j, b'-');
    }
    for (j, &op) in left.alignment.iter().enumerate() {
        if op == EDLIB_EDOP_INSERT {
            al.set(1, j + ms, b'-');
        } else {
            t_idx += 1;
            al.set(1, j + ms, target[t_idx as usize]);
        }
    }
    for (j, &op) in left.alignment.iter().enumerate() {
        if op == EDLIB_EDOP_DELETE {
            al.set(0, j + ms, b'-');
        } else {
            q_idx += 1;
            al.set(0, j + ms, query[q_idx as usize]);
        }
    }
    for j in 0..gl {
        al.set(0, j + ms + ll, b'-');
        t_idx += 1;
        al.set(1, j + ms + ll, target[t_idx as usize]);
    }
    for (j, &op) in right.alignment.iter().enumerate() {
        if op == EDLIB_EDOP_INSERT {
            al.set(1, j + ms + ll + gl, b'-');
        } else {
            t_idx += 1;
            al.set(1, j + ms + ll + gl, target[t_idx as usize]);
        }
    }
    for (j, &op) in right.alignment.iter().enumerate() {
        if op == EDLIB_EDOP_DELETE {
            al.set(0, j + ms + ll + gl, b'-');
        } else {
            q_idx += 1;
            al.set(0, j + ms + ll + gl, query[q_idx as usize]);
        }
    }
    for j in 0..me {
        t_idx += 1;
        al.set(1, j + ms + ll + gl + rl, target[t_idx as usize]);
        al.set(0, j + ms + ll + gl + rl, b'-');
    }
    al
}

/// splitAlign (insertions)
fn split_align(cons: &[u8], sv_ref: &[u8]) -> Option<Align> {
    let prefix = substr(sv_ref, 0, sv_ref.len() / 3);
    let cig_prefix = edlib::align(prefix, cons, -1, Mode::HW, Task::Path);
    let cs_start = edlib::infix_start(&cig_prefix);
    let suffix = substr_from(sv_ref, 2 * sv_ref.len() / 3);
    let cig_suffix = edlib::align(suffix, cons, -1, Mode::HW, Task::Path);
    let cs_end = edlib::infix_end(&cig_suffix);
    if cs_start >= cs_end {
        return None;
    }
    let cs = substr(cons, cs_start as usize, (cs_end - cs_start) as usize).to_vec();

    let cigar = edlib::align(sv_ref, &cs, -1, Mode::SHW, Task::Path);
    let dist_fwd = edit_distance_vec(sv_ref, &cs, &cigar);
    let mut sv_ref_rev = sv_ref.to_vec();
    let mut cs_rev = cs.clone();
    reverse_complement(&mut sv_ref_rev);
    reverse_complement(&mut cs_rev);
    let cigar_rev = edlib::align(&sv_ref_rev, &cs_rev, -1, Mode::SHW, Task::Path);
    let dist_rev = edit_distance_vec(&sv_ref_rev, &cs_rev, &cigar_rev);

    let sz = dist_fwd.len();
    let mut best_join: usize = 0;
    let mut i: usize = 1;
    while i < sz.wrapping_sub(1) {
        if dist_fwd[i].wrapping_add(dist_rev[sz - i - 2]) < dist_fwd[best_join].wrapping_add(dist_rev[sz - best_join - 2]) {
            best_join = i;
        }
        i += 1;
    }
    let left = substr(sv_ref, 0, best_join + 1);
    let right = substr_from(sv_ref, best_join + 1);
    let cig_left = edlib::align(left, cons, -1, Mode::HW, Task::Path);
    let left_end = edlib::infix_end(&cig_left);
    let cig_right = edlib::align(right, cons, -1, Mode::HW, Task::Path);
    let right_start = edlib::infix_start(&cig_right);
    if left_end.wrapping_add(15) >= right_start {
        return None;
    }
    Some(glue_alignment(sv_ref, cons, right_start - left_end - 1, &cig_left, &cig_right))
}

/// _consRefAlignment
pub fn cons_ref_alignment(cons: &[u8], sv_ref: &[u8], svt: i32) -> Option<Align> {
    if svt == 4 {
        let mut al = split_align(cons, sv_ref)?;
        for j in 0..al.cols {
            let tmp = al.get(0, j);
            al.set(0, j, al.get(1, j));
            al.set(1, j, tmp);
        }
        Some(al)
    } else {
        long_needle(cons, sv_ref, AlignConfig { horizontal: true, vertical: false }, &DnaScore::new(1, -1, -1, -1))
    }
}

/// _alignConsensus (realign = false)
fn align_consensus_inner(c: &Config, consensus: &[u8], sv_ref: &[u8], sv: &mut SVRecord, bp: &Breakpoint) -> bool {
    let al = match cons_ref_alignment(consensus, sv_ref, sv.svt) {
        Some(a) => a,
        None => return false,
    };
    let mut ad = AlignDescriptor::default();
    if !find_split(c, consensus, sv_ref, &al, &mut ad, sv.svt) {
        return false;
    }
    let (fgs, fge) = match coord_transform(c, sv_ref.len(), bp, &ad, sv.svt) {
        Some(v) => v,
        None => return false,
    };
    if is_translocation(sv.svt) || fgs < fge {
        if sv.sv_end - sv.sv_start <= c.indelsize && (sv.svt == 2 || sv.svt == 4) {
            let mut ref_vcf: Vec<u8> = Vec::new();
            let mut alt_vcf: Vec<u8> = Vec::new();
            let mut cpos: i32 = 0;
            let mut in_sv = false;
            for j in 0..al.cols {
                let a0 = al.get(0, j);
                let a1 = al.get(1, j);
                if a0 != b'-' {
                    cpos += 1;
                    if cpos == ad.c_start {
                        in_sv = true;
                    } else if cpos == ad.c_end {
                        in_sv = false;
                    }
                }
                if in_sv {
                    if a0 != b'-' {
                        alt_vcf.push(a0);
                    }
                    if a1 != b'-' {
                        ref_vcf.push(a1);
                    }
                }
            }
            let mut a = String::from_utf8_lossy(&ref_vcf).into_owned();
            a.push(',');
            a.push_str(&String::from_utf8_lossy(&alt_vcf));
            sv.alleles = a;
        }
        sv.precise = true;
        sv.sv_start = fgs as i32;
        sv.sv_end = fge as i32;
        sv.sr_align_quality = ad.perc_id;
        sv.ins_len = ad.c_end - ad.c_start - 1;
        sv.cons_bp = ad.c_start;
        sv.hom_len = 0.max(ad.hom_left + ad.hom_right - 2);
        let ci = ad.hom_left.max(ad.hom_right);
        sv.ciposlow = -ci;
        sv.ciposhigh = ci;
        sv.ciendlow = -ci;
        sv.ciendhigh = ci;
        true
    } else {
        false
    }
}

/// alignConsensus(c, hdr, seq, sndSeq, sv)
pub fn align_consensus(c: &Config, target_len: &[u32], seq: &RefSeq, snd_seq: Option<&RefSeq>, sv: &mut SVRecord) -> bool {
    if (sv.consensus.len() as i32) < (2 * c.minimum_flank_size + sv.ins_len) {
        return false;
    }
    let mut bp = Breakpoint::from_sv(sv);
    if sv.svt == 4 {
        let buffer = ((sv.consensus.len() as u64).wrapping_sub(sv.ins_len as i64 as u64) / 3) as i32;
        init_breakpoint(target_len, &mut bp, buffer.max(c.minimum_flank_size), sv.svt);
    } else {
        init_breakpoint(target_len, &mut bp, sv.consensus.len() as i32, sv.svt);
    }
    if bp.chr != bp.chr2 {
        bp.part1 = get_sv_ref(c, snd_seq.expect("second reference"), &bp, bp.chr2, sv.svt);
    }
    let sv_ref = get_sv_ref(c, seq, &bp, bp.chr, sv.svt);
    let consensus = sv.consensus.clone().into_bytes();
    align_consensus_inner(c, &consensus, &sv_ref, sv, &bp)
}
