//! Port of junction.h (split-read junction extraction, short-read path).

use crate::config::Config;
use crate::hts::*;
use crate::tags::*;
use crate::util::read_length;
use std::collections::HashMap;

pub type ReadBp = HashMap<u64, Vec<Junction>>;

/// _insertJunction
#[inline]
pub fn insert_junction(read_bp: &mut ReadBp, seed: u64, rec: &Record, rp: i32, sp: i32, scleft: bool) {
    let fw = rec.flag() & BAM_FREVERSE == 0;
    let mut read_start = rec.pos() as i32;
    if rec.flag() & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
        read_start = -1;
    }
    let seqlen = read_length(rec);
    if sp <= seqlen {
        let seqpos = if rec.flag() & BAM_FREVERSE != 0 { seqlen - sp } else { sp };
        read_bp.entry(seed).or_default().push(Junction { forward: fw, scleft, refidx: rec.tid(), rstart: read_start, refpos: rp, seqpos, qual: rec.qual() as u16 });
    }
}

fn select_read_start(v: &[Junction]) -> i32 {
    v.iter().find(|j| j.rstart != -1).map(|j| j.rstart).unwrap_or(-1)
}

fn select_primary_chr(v: &[Junction]) -> i32 {
    v.iter().find(|j| j.rstart != -1).map(|j| j.refidx).unwrap_or(-1)
}

#[inline]
fn avg_q(a: &Junction, b: &Junction) -> i32 {
    (a.qual as i32 + b.qual as i32) / 2
}

pub fn select_deletions(c: &Config, read_bp: &ReadBp, br: &mut [Vec<SRBamRecord>]) {
    for (&id, v) in read_bp {
        if v.len() <= 1 {
            continue;
        }
        let rst = select_read_start(v);
        for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if (v[j].seqpos - v[i].seqpos) as u32 > c.max_read_sep {
                    break;
                }
                if v[j].refidx == v[i].refidx && v[j].forward == v[i].forward && v[i].scleft != v[j].scleft {
                    let dellen = if v[i].forward {
                        if !v[i].scleft {
                            if v[i].refpos <= v[j].refpos {
                                (v[j].refpos - v[i].refpos) - (v[j].seqpos - v[i].seqpos)
                            } else {
                                0
                            }
                        } else if v[i].refpos <= v[j].refpos {
                            0
                        } else {
                            (v[i].refpos - v[j].refpos) + (v[j].seqpos - v[i].seqpos)
                        }
                    } else if v[i].scleft {
                        if v[i].refpos <= v[j].refpos {
                            0
                        } else {
                            (v[i].refpos - v[j].refpos) - (v[j].seqpos - v[i].seqpos)
                        }
                    } else if v[i].refpos <= v[j].refpos {
                        (v[j].refpos - v[i].refpos) + (v[j].seqpos - v[i].seqpos)
                    } else {
                        0
                    };
                    if dellen > c.min_ref_sep as i32 {
                        let q = avg_q(&v[i], &v[j]);
                        let sst = v[j].seqpos.min(v[i].seqpos);
                        let il = (v[j].seqpos - v[i].seqpos).abs();
                        if v[i].refpos <= v[j].refpos {
                            if !v[i].scleft && v[j].scleft {
                                br[2].push(SRBamRecord::new(v[i].refidx, v[i].refpos, v[j].refidx, v[j].refpos, rst, sst, q, il, id));
                            }
                        } else if v[i].scleft && !v[j].scleft {
                            br[2].push(SRBamRecord::new(v[j].refidx, v[j].refpos, v[i].refidx, v[i].refpos, rst, sst, q, il, id));
                        }
                        break;
                    }
                }
            }
        }
    }
}

pub fn select_duplications(c: &Config, read_bp: &ReadBp, br: &mut [Vec<SRBamRecord>]) {
    for (&id, v) in read_bp {
        if v.len() <= 1 {
            continue;
        }
        let rst = select_read_start(v);
        for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if (v[j].seqpos - v[i].seqpos) as u32 > c.max_read_sep {
                    break;
                }
                if v[j].refidx == v[i].refidx && v[j].forward == v[i].forward && v[i].scleft != v[j].scleft && (v[j].refpos - v[i].refpos).unsigned_abs() > c.min_ref_sep {
                    let q = avg_q(&v[i], &v[j]);
                    let sst = v[j].seqpos.min(v[i].seqpos);
                    let il = (v[j].seqpos - v[i].seqpos).abs();
                    if v[i].refpos <= v[j].refpos {
                        if v[i].scleft && !v[j].scleft {
                            br[3].push(SRBamRecord::new(v[i].refidx, v[i].refpos, v[j].refidx, v[j].refpos, rst, sst, q, il, id));
                        }
                    } else if !v[i].scleft && v[j].scleft {
                        br[3].push(SRBamRecord::new(v[j].refidx, v[j].refpos, v[i].refidx, v[i].refpos, rst, sst, q, il, id));
                    }
                }
            }
        }
    }
}

pub fn select_inversions(c: &Config, read_bp: &ReadBp, br: &mut [Vec<SRBamRecord>]) {
    for (&id, v) in read_bp {
        if v.len() <= 1 {
            continue;
        }
        let rst = select_read_start(v);
        for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if (v[j].seqpos - v[i].seqpos) as u32 > c.max_read_sep {
                    break;
                }
                if v[j].refidx == v[i].refidx && v[j].forward != v[i].forward && v[i].scleft == v[j].scleft && (v[j].refpos - v[i].refpos).unsigned_abs() > c.min_ref_sep {
                    let q = avg_q(&v[i], &v[j]);
                    let sst = v[j].seqpos.min(v[i].seqpos);
                    let il = (v[j].seqpos - v[i].seqpos).abs();
                    let t = if v[i].scleft { 1 } else { 0 };
                    if v[i].refpos <= v[j].refpos {
                        br[t].push(SRBamRecord::new(v[i].refidx, v[i].refpos, v[j].refidx, v[j].refpos, rst, sst, q, il, id));
                    } else {
                        br[t].push(SRBamRecord::new(v[j].refidx, v[j].refpos, v[i].refidx, v[i].refpos, rst, sst, q, il, id));
                    }
                }
            }
        }
    }
}

pub fn select_insertions(c: &Config, read_bp: &ReadBp, br: &mut [Vec<SRBamRecord>]) {
    for (&id, v) in read_bp {
        if v.len() <= 1 {
            continue;
        }
        let rst = select_read_start(v);
        'outer: for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if v[j].refidx == v[i].refidx && v[j].forward == v[i].forward && v[i].scleft != v[j].scleft && (v[j].refpos - v[i].refpos).unsigned_abs() < c.max_read_sep {
                    let isizelen = if v[i].forward {
                        if !v[i].scleft {
                            if v[i].refpos <= v[j].refpos {
                                (v[j].seqpos - v[i].seqpos) - (v[j].refpos - v[i].refpos)
                            } else {
                                (v[j].seqpos - v[i].seqpos) + (v[i].refpos - v[j].refpos)
                            }
                        } else {
                            0
                        }
                    } else if v[i].scleft {
                        if v[i].refpos <= v[j].refpos {
                            (v[j].seqpos - v[i].seqpos) + (v[j].refpos - v[i].refpos)
                        } else {
                            (v[j].seqpos - v[i].seqpos) - (v[i].refpos - v[j].refpos)
                        }
                    } else {
                        0
                    };
                    if isizelen > c.min_ref_sep as i32 && isizelen <= v[i].seqpos.max(v[j].seqpos) {
                        let q = avg_q(&v[i], &v[j]);
                        let sst = v[j].seqpos.min(v[i].seqpos);
                        if v[i].refpos <= v[j].refpos {
                            br[4].push(SRBamRecord::new(v[i].refidx, v[i].refpos, v[j].refidx, v[j].refpos, rst, sst, q, isizelen, id));
                        } else {
                            br[4].push(SRBamRecord::new(v[j].refidx, v[j].refpos, v[i].refidx, v[i].refpos, rst, sst, q, isizelen, id));
                        }
                        // `break` only leaves the inner j-loop
                        continue 'outer;
                    }
                }
            }
        }
    }
}

pub fn select_translocations(c: &Config, read_bp: &ReadBp, br: &mut [Vec<SRBamRecord>]) {
    for (&id, v) in read_bp {
        if v.len() <= 1 {
            continue;
        }
        let rst = select_read_start(v);
        let primary = select_primary_chr(v);
        for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if (v[j].seqpos - v[i].seqpos) as u32 > c.max_read_sep {
                    break;
                }
                if v[j].refidx != v[i].refidx {
                    let (c1, c2) = if v[i].refidx < v[j].refidx { (i, j) } else { (j, i) };
                    let q = avg_q(&v[i], &v[j]);
                    let mk = || {
                        let mut r = SRBamRecord::new(v[c2].refidx, v[c2].refpos, v[c1].refidx, v[c1].refpos, rst, v[j].seqpos.min(v[i].seqpos), q, (v[j].seqpos - v[i].seqpos).abs(), id);
                        r.primary_chr = primary;
                        r
                    };
                    if v[c1].forward == v[c2].forward {
                        if v[c1].scleft != v[c2].scleft {
                            if v[c1].scleft {
                                br[(DELLY_SVT_TRANS + 2) as usize].push(mk());
                            } else {
                                br[(DELLY_SVT_TRANS + 3) as usize].push(mk());
                            }
                        }
                    } else if v[c1].scleft == v[c2].scleft {
                        if v[c1].scleft {
                            br[(DELLY_SVT_TRANS + 1) as usize].push(mk());
                        } else {
                            br[DELLY_SVT_TRANS as usize].push(mk());
                        }
                    }
                }
            }
        }
    }
}
