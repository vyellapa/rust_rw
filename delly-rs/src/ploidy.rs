//! Port of ploidy.h (sex chromosome handling for `delly sr`).

use crate::hts::*;

#[derive(Clone, Debug)]
pub struct SexModel {
    pub x_tid: i32,
    pub y_tid: i32,
    pub par_x: Vec<(i32, i32)>,
    pub par_y: Vec<(i32, i32)>,
    pub sex: Vec<u8>,
}

impl Default for SexModel {
    fn default() -> Self {
        SexModel { x_tid: -1, y_tid: -1, par_x: Vec::new(), par_y: Vec::new(), sex: Vec::new() }
    }
}

fn in_regions(reg: &[(i32, i32)], pos: i32) -> bool {
    reg.iter().any(|&(a, b)| pos >= a && pos < b)
}

pub fn ploidy(sm: &SexModel, sex: u8, tid: i32, pos: i32) -> u8 {
    if sex == 0 {
        return 2;
    }
    if tid == sm.x_tid {
        if sex == 1 && !in_regions(&sm.par_x, pos) {
            return 1;
        }
        return 2;
    }
    if tid == sm.y_tid {
        if sex == 2 {
            return 0;
        }
        if in_regions(&sm.par_y, pos) {
            return 2;
        }
        return 1;
    }
    2
}

pub fn sv_ploidy(sm: &SexModel, sex: u8, tid: i32, pos: i32, tid2: i32, pos2: i32) -> u8 {
    let p1 = ploidy(sm, sex, tid, pos);
    let p2 = ploidy(sm, sex, tid2, pos2);
    if p1 == 0 || p2 == 0 {
        return 0;
    }
    p1.max(p2)
}

fn par_regions(xlen: u32, sm: &mut SexModel) -> bool {
    sm.par_x.clear();
    sm.par_y.clear();
    match xlen {
        156040895 => {
            sm.par_x = vec![(10000, 2781479), (155701382, 156030895)];
            sm.par_y = vec![(10000, 2781479), (56887902, 57217415)];
        }
        155270560 => {
            sm.par_x = vec![(60000, 2699520), (154931043, 155260560)];
            sm.par_y = vec![(10000, 2649520), (59034049, 59363566)];
        }
        154259566 => {
            sm.par_x = vec![(0, 2394410), (153925834, 154259566)];
            sm.par_y = vec![(0, 2458320), (62122809, 62460029)];
        }
        _ => return false,
    }
    true
}

/// _sexChromosomes (BAM header)
pub fn sex_chromosomes(hdr: &Header, sm: &mut SexModel) -> bool {
    sm.x_tid = -1;
    sm.y_tid = -1;
    for (i, n) in hdr.names.iter().enumerate() {
        if n == "chrX" || n == "X" {
            sm.x_tid = i as i32;
        } else if n == "chrY" || n == "Y" {
            sm.y_tid = i as i32;
        }
    }
    if sm.x_tid == -1 {
        return false;
    }
    if !par_regions(hdr.lens[sm.x_tid as usize], sm) {
        sm.x_tid = -1;
        sm.y_tid = -1;
        return false;
    }
    true
}

/// _inferSexFromReads; returns (sex, xRatio, yRatio)
pub fn infer_sex_from_reads(path: &str, genome: &str, sm: &SexModel) -> (u8, f64, f64) {
    let mut x_ratio = -1.0;
    let mut y_ratio = -1.0;
    if sm.x_tid == -1 {
        return (0, x_ratio, y_ratio);
    }
    let mut bam = match IndexedBam::open(path, Some(genome)) {
        Some(b) => b,
        None => return (0, x_ratio, y_ratio),
    };
    bam.set_required_fields(SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ);
    let hdr = bam.header.clone();
    let nwin: i32 = 20;
    let wlen: i32 = 200000;
    let mut auto_iv: Vec<(i32, (i32, i32))> = Vec::new();
    for (ri, name) in hdr.names.iter().enumerate() {
        let mut t: &str = name;
        if t.len() >= 3 && &t[..3] == "chr" {
            t = &t[3..];
        }
        let is_number = !t.is_empty() && t.bytes().all(|c| c.is_ascii_digit());
        if is_number && hdr.lens[ri] > (2 * nwin * wlen) as u32 {
            auto_iv.push((ri as i32, (wlen, hdr.lens[ri] as i32 - wlen)));
        }
    }
    if auto_iv.is_empty() {
        return (0, x_ratio, y_ratio);
    }
    let mut x_beg = wlen;
    let mut x_end = hdr.lens[sm.x_tid as usize] as i32 - wlen;
    if !sm.par_x.is_empty() {
        x_beg = sm.par_x[0].1 + wlen;
        x_end = sm.par_x[sm.par_x.len() - 1].0 - wlen;
    }
    let mut y_beg = -1;
    let mut y_end = -1;
    if sm.y_tid != -1 {
        y_beg = if sm.par_y.is_empty() { wlen } else { sm.par_y[0].1 + wlen };
        y_end = (0.45 * hdr.lens[sm.y_tid as usize] as f64) as i32;
    }
    let mut cnt: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for cls in 0..3 {
        let iv: Vec<(i32, (i32, i32))> = if cls == 0 {
            auto_iv.clone()
        } else if cls == 1 {
            vec![(sm.x_tid, (x_beg, x_end))]
        } else {
            if sm.y_tid == -1 {
                continue;
            }
            vec![(sm.y_tid, (y_beg, y_end))]
        };
        for &(tid, (b, e)) in &iv {
            let span = e - b;
            if span < nwin * wlen {
                continue;
            }
            let step = span / nwin;
            for w in 0..nwin {
                let beg = b + w * step;
                let mut c: u32 = 0;
                let ok = bam.for_each_in(tid, beg as i64, (beg + wlen) as i64, |rec| {
                    if rec.flag() & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP) != 0 {
                        return true;
                    }
                    if rec.qual() < 20 {
                        return true;
                    }
                    c += 1;
                    true
                });
                if !ok {
                    continue;
                }
                cnt[cls].push(c as f64);
            }
        }
    }
    if cnt[0].len() < 10 || cnt[1].len() < 5 {
        return (0, x_ratio, y_ratio);
    }
    let mut med = [0.0f64; 3];
    for cls in 0..3 {
        if cnt[cls].is_empty() {
            continue;
        }
        cnt[cls].sort_by(|a, b| a.partial_cmp(b).unwrap());
        med[cls] = cnt[cls][cnt[cls].len() / 2];
    }
    if med[0] <= 0.0 {
        return (0, x_ratio, y_ratio);
    }
    x_ratio = med[1] / med[0];
    y_ratio = if cnt[2].is_empty() { -1.0 } else { med[2] / med[0] };
    if x_ratio < 0.75 {
        if y_ratio < 0.0 || y_ratio > 0.1 {
            return (1, x_ratio, y_ratio);
        }
    } else if y_ratio < 0.0 || y_ratio < 0.15 {
        return (2, x_ratio, y_ratio);
    }
    (0, x_ratio, y_ratio)
}

/// _parseSex; returns Some(autoSex) or None on error
pub fn parse_sex(sample_names: &[String], sex_arg: &str, sm: &mut SexModel) -> Option<bool> {
    let n = sample_names.len();
    sm.sex = vec![0; n];
    match sex_arg {
        "auto" => return Some(true),
        "none" => return Some(false),
        "male" => sm.sex = vec![1; n],
        "female" => sm.sex = vec![2; n],
        _ => {
            let content = match std::fs::read(sex_arg) {
                Ok(c) => c,
                Err(_) => {
                    eprintln!("Sex file cannot be opened: {}", sex_arg);
                    return None;
                }
            };
            let text = String::from_utf8_lossy(&content);
            let mut map: std::collections::BTreeMap<String, u8> = std::collections::BTreeMap::new();
            for line in text.split('\n') {
                // boost::split with token_compress_on: leading delimiter still yields an empty first token
                let mut tokens: Vec<&str> = Vec::new();
                let mut cur_start = 0usize;
                let bytes = line.as_bytes();
                let mut i = 0;
                while i <= bytes.len() {
                    if i == bytes.len() || matches!(bytes[i], b'\t' | b' ' | b',') {
                        tokens.push(&line[cur_start..i]);
                        let mut j = i;
                        while j < bytes.len() && matches!(bytes[j], b'\t' | b' ' | b',') {
                            j += 1;
                        }
                        if i == bytes.len() {
                            break;
                        }
                        cur_start = j;
                        i = j;
                        if j == bytes.len() {
                            tokens.push("");
                            break;
                        }
                        continue;
                    }
                    i += 1;
                }
                if tokens.len() < 2 {
                    continue;
                }
                let s = tokens[1].to_ascii_lowercase();
                if s == "male" || s == "m" || s == "1" {
                    map.insert(tokens[0].to_string(), 1);
                } else if s == "female" || s == "f" || s == "2" {
                    map.insert(tokens[0].to_string(), 2);
                }
            }
            for (i, name) in sample_names.iter().enumerate() {
                if let Some(&v) = map.get(name) {
                    sm.sex[i] = v;
                }
            }
        }
    }
    Some(false)
}

pub fn sex_name(sex: u8) -> &'static str {
    match sex {
        1 => "male",
        2 => "female",
        _ => "unknown",
    }
}
