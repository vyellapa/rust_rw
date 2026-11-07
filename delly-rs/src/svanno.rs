//! Port of svanno.h: breakpoint homology, MEI/NUMT and tandem-repeat annotation.

use crate::config::Config;
use crate::edlib::{self, Mode, Task, EDLIB_STATUS_OK};
use crate::hts::RefSeq;
use crate::mei;
use crate::tags::SVRecord;
use crate::util::{reverse_complement, to_upper};

/// detectTandemRepeat (autocorrelation)
pub fn detect_tandem_repeat(seq: &[u8], max_period: i32, min_fraction: f32) -> (i32, f32) {
    let n = seq.len() as i32;
    if n < 10 {
        return (0, 0.0);
    }
    let pmax = max_period.min(n / 2);
    for p in 1..=pmax {
        let mut matches: i32 = 0;
        for i in p..n {
            if seq[i as usize] == seq[(i - p) as usize] {
                matches += 1;
            }
        }
        if matches as f32 / (n - p) as f32 >= min_fraction {
            return (p, n as f32 / p as f32);
        }
    }
    (0, 0.0)
}

#[inline]
fn upc(c: u8) -> u8 {
    c.to_ascii_uppercase()
}

/// annotateSV
pub fn annotate_sv(c: &Config, target_len: &[u32], seq: &RefSeq, sv: &mut SVRecord) {
    let chr_len = target_len[sv.chr as usize] as i32;
    let max_edits = 2;
    let at = |i: i32| -> u8 { seq.at(i as usize) };

    let mut alt_str: Vec<u8> = Vec::new();
    let mut ins_seq: Vec<u8> = Vec::new();
    if sv.svt == 4 && !sv.alleles.is_empty() {
        if let Some(comma) = sv.alleles.find(',') {
            let alt = sv.alleles.as_bytes()[comma + 1..].to_vec();
            if alt.len() > 1 && alt[0] != b'<' {
                alt_str = to_upper(&alt);
                ins_seq = alt_str[1..].to_vec();
                if ins_seq.len() > 15000 {
                    ins_seq.truncate(15000);
                }
            }
        }
    }

    if sv.svt == 4 {
        if !alt_str.is_empty() {
            let limit = 100.min(alt_str.len() as i32).min(chr_len - (sv.sv_start - 1));
            let mut edits = 0;
            for h in 1..=limit {
                if upc(alt_str[(h - 1) as usize]) != upc(at(sv.sv_start - 1 + h - 1)) {
                    edits += 1;
                    if edits > max_edits {
                        break;
                    }
                }
                sv.anno.hom_len = h;
            }
        }
    } else if sv.sv_end > sv.sv_start {
        let max_h_del = 10000;
        let min_identity = 0.90f64;
        let stop_identity = 0.75f64;
        let early_stop_gap = 100;
        let half_sv_len = (sv.sv_end - sv.sv_start) / 2;
        let bwd_limit = max_h_del.min(half_sv_len).min(sv.sv_start);
        if bwd_limit > 0 {
            let mut edits = 0;
            let mut last_good_h = 0;
            for h in 1..=bwd_limit {
                if upc(at(sv.sv_start - h)) != upc(at(sv.sv_end - h)) {
                    edits += 1;
                }
                let id = 1.0 - edits as f64 / h as f64;
                if id >= min_identity {
                    sv.anno.hom_len = h;
                    last_good_h = h;
                } else if h - last_good_h > early_stop_gap && id < stop_identity {
                    break;
                }
            }
        }
        let fwd_limit = 0.max(max_h_del.min(half_sv_len).min(chr_len - sv.sv_end - 1));
        if fwd_limit > 0 {
            let mut edits = 0;
            let mut last_good_h = 0;
            for h in 1..=fwd_limit {
                if upc(at(sv.sv_start + h)) != upc(at(sv.sv_end + h)) {
                    edits += 1;
                }
                let id = 1.0 - edits as f64 / h as f64;
                if id >= min_identity {
                    if h > sv.anno.hom_len {
                        sv.anno.hom_len = h;
                    }
                    last_good_h = h;
                } else if h - last_good_h > early_stop_gap && id < stop_identity {
                    break;
                }
            }
        }
    }

    if !ins_seq.is_empty() {
        let mut best_id = c.mei_min_frac as f64;
        let mut best_type = 0;
        let mut best_fwd = 0.0f64;
        let mut best_rev = 0.0f64;
        const MIN_INS_LEN: [i32; 7] = [0, 100, 150, 1000, 100, 300, 300];
        for seqtype in 1..=6usize {
            if (ins_seq.len() as i32) < MIN_INS_LEN[seqtype] {
                continue;
            }
            let mut templ: Vec<u8> = match seqtype {
                1 => mei::ALU,
                2 => mei::LINE1,
                3 => mei::SVA,
                4 => mei::NUMT,
                5 => mei::SOLOLTR,
                _ => mei::HERVK,
            }
            .as_bytes()
            .to_vec();
            if seqtype < 4 {
                templ.extend_from_slice(mei::POLYA.as_bytes());
            }
            let mut templ_rc = templ.clone();
            reverse_complement(&mut templ_rc);
            let long_ins = ins_seq.len() > templ.len();
            let (fq, ft) = if long_ins { (&templ, &ins_seq) } else { (&ins_seq, &templ) };
            let mut fwd_id = 0.0f64;
            let r = edlib::align(fq, ft, -1, Mode::HW, Task::Distance);
            if r.status == EDLIB_STATUS_OK && r.edit_distance >= 0 {
                fwd_id = 1.0 - r.edit_distance as f64 / fq.len() as i32 as f64;
            }
            let (rq, rt) = if long_ins { (&templ_rc, &ins_seq) } else { (&ins_seq, &templ_rc) };
            let mut rev_id = 0.0f64;
            let r = edlib::align(rq, rt, -1, Mode::HW, Task::Distance);
            if r.status == EDLIB_STATUS_OK && r.edit_distance >= 0 {
                rev_id = 1.0 - r.edit_distance as f64 / rq.len() as i32 as f64;
            }
            let coverage = if long_ins { fq.len() as i32 as f64 / ft.len() as i32 as f64 } else { 1.0 };
            let effective = fwd_id.max(rev_id).min(coverage);
            if effective > best_id {
                best_id = effective;
                best_type = seqtype as i32;
                best_fwd = fwd_id;
                best_rev = rev_id;
            }
        }
        if best_type > 0 {
            sv.anno.seq_type = best_type;
            sv.anno.is_rc = best_rev > best_fwd;
        } else if ins_seq.len() >= 10 {
            let mut period = 0;
            let mut copies = 0.0f32;
            let (p1, c1) = detect_tandem_repeat(&ins_seq, 100, c.tr_min_frac);
            if p1 > 0 {
                period = p1;
                copies = c1;
            } else {
                let flank = 200;
                let mut side = 0;
                while side < 2 && period == 0 {
                    let ref_win: Vec<u8> = if side == 0 {
                        let rbeg = 0.max(sv.sv_start - flank);
                        if sv.sv_start - rbeg < 40 {
                            side += 1;
                            continue;
                        }
                        to_upper(&seq.seq()[rbeg as usize..sv.sv_start as usize])
                    } else {
                        let rend = chr_len.min(sv.sv_start + flank);
                        if rend - sv.sv_start < 40 {
                            side += 1;
                            continue;
                        }
                        to_upper(&seq.seq()[sv.sv_start as usize..rend as usize])
                    };
                    let (pr, _cr) = detect_tandem_repeat(&ref_win, 100, c.tr_min_frac);
                    if pr <= 0 || (ref_win.len() as i32) < pr {
                        side += 1;
                        continue;
                    }
                    let unit: Vec<u8> = if side == 0 { ref_win[ref_win.len() - pr as usize..].to_vec() } else { ref_win[..pr as usize].to_vec() };
                    let mut templ: Vec<u8> = Vec::with_capacity(ins_seq.len() + pr as usize);
                    while (templ.len() as i32) < ins_seq.len() as i32 + pr {
                        templ.extend_from_slice(&unit);
                    }
                    let r = edlib::align(&ins_seq, &templ, -1, Mode::HW, Task::Distance);
                    let id = if r.status == EDLIB_STATUS_OK && r.edit_distance >= 0 { 1.0 - r.edit_distance as f64 / ins_seq.len() as f64 } else { 0.0 };
                    if id >= 0.70 {
                        period = pr;
                        copies = ins_seq.len() as f32 / pr as f32;
                    }
                    side += 1;
                }
            }
            if period > 0 {
                sv.anno.seq_type = 7;
                sv.anno.tr_period = period;
                sv.anno.tr_copies = copies;
            }
        }
    } else if sv.svt == 2 {
        let del_len = sv.sv_end - sv.sv_start;
        if (10..=50000).contains(&del_len) {
            let del_seq = to_upper(&seq.seq()[sv.sv_start as usize..sv.sv_end as usize]);
            let (period, copies) = detect_tandem_repeat(&del_seq, 100, c.tr_min_frac);
            if period > 0 {
                sv.anno.seq_type = 7;
                sv.anno.tr_period = period;
                sv.anno.tr_copies = copies;
            }
        }
    }
}
