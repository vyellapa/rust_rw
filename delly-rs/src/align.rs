//! Port of align.h / gotoh.h / needle.h (the parts used by `delly sr`).

use crate::config::DnaScore;
use crate::util::reverse_complement;

/// 2-D character alignment (boost::multi_array<char, 2>), row-major
#[derive(Clone, Debug, Default)]
pub struct Align {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<u8>,
}

impl Align {
    pub fn new(rows: usize, cols: usize) -> Self {
        Align { rows, cols, data: vec![0u8; rows * cols] }
    }
    pub fn from_seq(s: &[u8]) -> Self {
        Align { rows: 1, cols: s.len(), data: s.to_vec() }
    }
    #[inline]
    pub fn get(&self, r: usize, c: usize) -> u8 {
        self.data[r * self.cols + c]
    }
    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: u8) {
        self.data[r * self.cols + c] = v;
    }
}

/// AlignConfig<horizontal, vertical>: free end gaps
#[derive(Clone, Copy)]
pub struct AlignConfig {
    pub horizontal: bool,
    pub vertical: bool,
}

#[inline]
fn vertical_gap(ac: AlignConfig, i: usize, iend: usize, cost: i32) -> i32 {
    if ac.vertical && (i == 0 || i == iend) {
        0
    } else {
        cost
    }
}

#[inline]
fn horizontal_gap(ac: AlignConfig, i: usize, iend: usize, cost: i32) -> i32 {
    if ac.horizontal && (i == 0 || i == iend) {
        0
    } else {
        cost
    }
}

/// Alignment profile (6 x cols), float as in gotoh.h
struct Profile {
    cols: usize,
    p: Vec<f32>,
}

impl Profile {
    #[inline]
    fn at(&self, k: usize, j: usize) -> f32 {
        self.p[k * self.cols + j]
    }
}

/// _createProfile(multi_array)
fn create_profile(a: &Align) -> Profile {
    let cols = a.cols;
    let mut p = vec![0f32; 6 * cols];
    let mut first = vec![-1i32; a.rows];
    let mut last = vec![cols as i32; a.rows];
    for i in 0..a.rows {
        for j in 0..cols {
            let ch = a.get(i, j);
            if first[i] == -1 && ch != b'-' {
                first[i] = j as i32;
            }
            if first[i] != -1 && ch != b'-' {
                last[i] = j as i32;
            }
        }
    }
    for j in 0..cols {
        let mut sum: i32 = 0;
        for i in 0..a.rows {
            if first[i] <= j as i32 && j as i32 <= last[i] {
                sum += 1;
                match a.get(i, j) {
                    b'A' | b'a' => p[j] += 1.0,
                    b'C' | b'c' => p[cols + j] += 1.0,
                    b'G' | b'g' => p[2 * cols + j] += 1.0,
                    b'T' | b't' => p[3 * cols + j] += 1.0,
                    b'N' | b'n' => p[4 * cols + j] += 1.0,
                    b'-' => p[5 * cols + j] += 1.0,
                    _ => sum -= 1,
                }
            }
        }
        for k in 0..6 {
            p[k * cols + j] /= sum as f32;
        }
    }
    Profile { cols, p }
}

/// _score for multi_arrays. clang contracts `score += p1*p2*w` into fmadd(p1*p2, w, score).
#[inline(always)]
fn profile_score(a1: &Align, a2: &Align, p: Option<&(Profile, Profile)>, row: usize, col: usize, sc: &DnaScore) -> i32 {
    if a1.rows == 1 && a2.rows == 1 {
        if a1.get(0, row) == a2.get(0, col) {
            sc.match_
        } else {
            sc.mismatch
        }
    } else {
        let (p1, p2) = p.unwrap();
        let mut score: f32 = 0.0;
        let wm = sc.match_ as f32;
        let wx = sc.mismatch as f32;
        for k1 in 0..5 {
            let a = p1.at(k1, row);
            for k2 in 0..5 {
                let w = if k1 == k2 { wm } else { wx };
                score = (a * p2.at(k2, col)).mul_add(w, score);
            }
        }
        score as i32
    }
}

/// gotoh on two multi_arrays producing a merged multi_array
pub fn gotoh(a1: &Align, a2: &Align, ac: AlignConfig, sc: &DnaScore) -> Align {
    // Dispatch once per alignment, rather than once per mul_add in every DP cell.
    // Both paths retain fused rounding; no fast-math or reassociation is enabled.
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::is_x86_feature_detected!("fma") {
        // SAFETY: the CPU/OS feature check guards the only entry to this path.
        return unsafe { gotoh_fma(a1, a2, ac, sc) };
    }
    gotoh_impl(a1, a2, ac, sc)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "fma")]
unsafe fn gotoh_fma(a1: &Align, a2: &Align, ac: AlignConfig, sc: &DnaScore) -> Align {
    gotoh_impl(a1, a2, ac, sc)
}

#[inline(always)]
fn gotoh_impl(a1: &Align, a2: &Align, ac: AlignConfig, sc: &DnaScore) -> Align {
    let m = a1.cols;
    let n = a2.cols;
    let mut s = vec![0i32; n + 1];
    let mut v = vec![0i32; n + 1];
    let mut newhoz: i32 = 0;
    let mut prevsub: i32 = 0;
    let mf = n + 1;
    let tot = (m + 1) * (n + 1);
    // bit1..bit4 packed into one byte per cell
    let mut bits = vec![0u8; tot];
    const B1: u8 = 1;
    const B2: u8 = 2;
    const B3: u8 = 4;
    const B4: u8 = 8;
    let prof = if a1.rows != 1 || a2.rows != 1 { Some((create_profile(a1), create_profile(a2))) } else { None };
    let inf = sc.inf;
    for row in 0..=m {
        for col in 0..=n {
            if row == 0 && col == 0 {
                s[0] = 0;
                v[0] = -inf;
                newhoz = -inf;
                bits[0] |= B1 | B2;
            } else if row == 0 {
                v[col] = -inf;
                s[col] = horizontal_gap(ac, 0, m, sc.go + col as i32 * sc.ge);
                newhoz = horizontal_gap(ac, 0, m, sc.go + col as i32 * sc.ge);
                bits[col] |= B3;
            } else if col == 0 {
                newhoz = -inf;
                s[0] = vertical_gap(ac, 0, n, sc.go + row as i32 * sc.ge);
                if row - 1 == 0 {
                    prevsub = 0;
                } else {
                    prevsub = vertical_gap(ac, 0, n, sc.go + (row as i32 - 1) * sc.ge);
                }
                v[0] = vertical_gap(ac, 0, n, sc.go + row as i32 * sc.ge);
                bits[row * mf] |= B4;
            } else {
                let prevhoz = newhoz;
                let prevver = v[col];
                let prevprevsub = prevsub;
                prevsub = s[col];
                newhoz = (s[col - 1] + horizontal_gap(ac, row, m, sc.go + sc.ge)).max(prevhoz + horizontal_gap(ac, row, m, sc.ge));
                v[col] = (prevsub + vertical_gap(ac, col, n, sc.go + sc.ge)).max(prevver + vertical_gap(ac, col, n, sc.ge));
                s[col] = (prevprevsub + profile_score(a1, a2, prof.as_ref(), row - 1, col - 1, sc)).max(newhoz).max(v[col]);
                let idx = row * mf + col;
                if s[col] == newhoz {
                    bits[idx] |= B3;
                } else if s[col] == v[col] {
                    bits[idx] |= B4;
                }
                if newhoz != prevhoz + horizontal_gap(ac, row, m, sc.ge) {
                    bits[idx] |= B1;
                }
                if v[col] != prevver + vertical_gap(ac, col, n, sc.ge) {
                    bits[idx] |= B2;
                }
            }
        }
    }
    // Trace-back
    let mut row = m;
    let mut col = n;
    let mut last = b's';
    let mut btr: Vec<u8> = Vec::new();
    while row > 0 || col > 0 {
        let b = bits[row * mf + col];
        if last == b's' {
            if b & B3 != 0 {
                last = b'h';
            } else if b & B4 != 0 {
                last = b'v';
            } else {
                row -= 1;
                col -= 1;
                btr.push(b's');
            }
        } else if last == b'h' {
            if b & B1 != 0 {
                last = b's';
            }
            col -= 1;
            btr.push(b'h');
        } else {
            if b & B2 != 0 {
                last = b's';
            }
            row -= 1;
            btr.push(b'v');
        }
    }
    create_alignment_ma(&btr, a1, a2)
}

/// _createAlignment for multi_arrays
fn create_alignment_ma(trace: &[u8], a1: &Align, a2: &Align) -> Align {
    let num_n = a1.rows;
    let num_m = a2.rows;
    let mut al = Align::new(num_n + num_m, trace.len());
    let mut row = 0usize;
    let mut col = 0usize;
    for (ai, &t) in trace.iter().rev().enumerate() {
        if t == b's' {
            for i in 0..num_n {
                al.set(i, ai, a1.get(i, row));
            }
            for i in 0..num_m {
                al.set(num_n + i, ai, a2.get(i, col));
            }
            row += 1;
            col += 1;
        } else if t == b'h' {
            for i in 0..num_n {
                al.set(i, ai, b'-');
            }
            for i in 0..num_m {
                al.set(num_n + i, ai, a2.get(i, col));
            }
            col += 1;
        } else {
            for i in 0..num_n {
                al.set(i, ai, a1.get(i, row));
            }
            for i in 0..num_m {
                al.set(num_n + i, ai, b'-');
            }
            row += 1;
        }
    }
    al
}

/// _createLocalAlignment(trace, s1, s2, align, 0, 0) for strings
fn create_alignment_str(trace: &[u8], s1: &[u8], s2: &[u8]) -> Align {
    let mut al = Align::new(2, trace.len());
    let mut row = 0usize;
    let mut col = 0usize;
    for (ai, &t) in trace.iter().rev().enumerate() {
        if t == b's' {
            al.set(0, ai, s1[row]);
            row += 1;
            al.set(1, ai, s2[col]);
            col += 1;
        } else if t == b'h' {
            al.set(0, ai, b'-');
            al.set(1, ai, s2[col]);
            col += 1;
        } else {
            al.set(0, ai, s1[row]);
            row += 1;
            al.set(1, ai, b'-');
        }
    }
    al
}

/// longestHomology (banded edit distance)
pub fn longest_homology(s1: &[u8], s2: &[u8], score_threshold: i32) -> i32 {
    let m = s1.len() as i32;
    let n = s2.len() as i32;
    let k = score_threshold.abs();
    let w = (2 * k + 1) as usize;
    let mut prev = vec![0i32; w];
    let mut cur = vec![0i32; w];
    for h in 0..=k {
        prev[(h + k) as usize] = -h;
    }
    for row in 1..=m {
        let mut best_col = score_threshold - 1;
        for h in -k..=k {
            let col = row + h;
            if col == 0 {
                cur[(h + k) as usize] = -row;
            } else if col >= 1 && col <= n {
                let mut val = prev[(h + k) as usize] + if s1[(row - 1) as usize] == s2[(col - 1) as usize] { 0 } else { -1 };
                if h + 1 <= k {
                    val = val.max(prev[(h + 1 + k) as usize] - 1);
                }
                if h - 1 >= -k {
                    val = val.max(cur[(h - 1 + k) as usize] - 1);
                }
                cur[(h + k) as usize] = val;
                if val > best_col {
                    best_col = val;
                }
            }
        }
        if best_col < score_threshold {
            return row - 1;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    0
}

const NEG: i32 = -(1 << 28);

#[derive(Default)]
struct SplitDPMatrix {
    lo: Vec<i32>,
    hi: Vec<i32>,
    off: Vec<usize>,
    val: Vec<i32>,
}

impl SplitDPMatrix {
    fn reset(&mut self, m: i32) {
        self.lo.clear();
        self.lo.resize((m + 1) as usize, 1);
        self.hi.clear();
        self.hi.resize((m + 1) as usize, 0);
        self.off.clear();
        self.off.resize((m + 1) as usize, 0);
        self.val.clear();
    }
    #[inline]
    fn empty(&self, r: i32) -> bool {
        self.lo[r as usize] > self.hi[r as usize]
    }
    #[inline]
    fn get(&self, r: i32, c: i32) -> i32 {
        let ru = r as usize;
        if c < self.lo[ru] || c > self.hi[ru] {
            return NEG;
        }
        self.val[self.off[ru] + (c - self.lo[ru]) as usize]
    }
    fn trim(&mut self, r: i32, first: i32, last: i32, thr: i32) {
        let ru = r as usize;
        let start = self.off[ru];
        let mut lokeep = first;
        while lokeep <= last && self.val[start + (lokeep - first) as usize] < thr {
            lokeep += 1;
        }
        let mut hikeep = last;
        while hikeep >= lokeep && self.val[start + (hikeep - first) as usize] < thr {
            hikeep -= 1;
        }
        if lokeep > hikeep {
            self.lo[ru] = 1;
            self.hi[ru] = 0;
            self.val.truncate(start);
        } else {
            self.lo[ru] = lokeep;
            self.hi[ru] = hikeep;
            self.off[ru] = start + (lokeep - first) as usize;
            self.val.resize(start + (hikeep - first) as usize + 1, 0);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn split_fill(s1: &[u8], s2: &[u8], ac: AlignConfig, sc: &DnaScore, lower_bound: i32, gain_max: i32, prune: bool, dp: &mut SplitDPMatrix) {
    let m = s1.len() as i32;
    let n = s2.len() as i32;
    dp.reset(m);
    dp.val.reserve(4 * (n as usize + 1));
    let mut thr = if prune { lower_bound - m * gain_max } else { NEG };
    dp.off[0] = dp.val.len();
    let mut v: i32 = 0;
    dp.val.push(v);
    for _col in 1..=n {
        v += horizontal_gap(ac, 0, m as usize, sc.ge);
        dp.val.push(v);
    }
    dp.trim(0, 0, n, thr);
    for row in 1..=m {
        if dp.empty(row - 1) {
            break;
        }
        thr = if prune { lower_bound - (m - row) * gain_max } else { NEG };
        let plo = dp.lo[(row - 1) as usize];
        let phi = dp.hi[(row - 1) as usize];
        let poff = dp.off[(row - 1) as usize];
        let hg = horizontal_gap(ac, row as usize, m as usize, sc.ge);
        let c1 = s1[(row - 1) as usize];
        dp.off[row as usize] = dp.val.len();
        let mut left = dp.val[poff] + vertical_gap(ac, plo as usize, n as usize, sc.ge);
        dp.val.push(left);
        for col in (plo + 1)..=phi {
            let diag = dp.val[poff + (col - 1 - plo) as usize] + if c1 == s2[(col - 1) as usize] { sc.match_ } else { sc.mismatch };
            let up = dp.val[poff + (col - plo) as usize] + vertical_gap(ac, col as usize, n as usize, sc.ge);
            left = diag.max(up).max(left + hg);
            dp.val.push(left);
        }
        let mut col = phi + 1;
        if col <= n {
            let diag = dp.val[poff + (phi - plo) as usize] + if c1 == s2[(col - 1) as usize] { sc.match_ } else { sc.mismatch };
            left = diag.max(left + hg);
            dp.val.push(left);
            col += 1;
            while col <= n && left + hg >= thr {
                left += hg;
                dp.val.push(left);
                col += 1;
            }
        }
        dp.trim(row, plo, col - 1, thr);
    }
}

#[inline]
fn is_acgtn(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'N')
}

#[allow(clippy::too_many_arguments)]
fn split_traceback(dp: &SplitDPMatrix, mut rr: i32, mut cc: i32, m: i32, n: i32, ac: AlignConfig, sc: &DnaScore, trace: &mut Vec<u8>) -> bool {
    while rr > 0 || cc > 0 {
        let v = dp.get(rr, cc);
        if v == NEG {
            return false;
        }
        if rr > 0 && v == dp.get(rr - 1, cc) + vertical_gap(ac, cc as usize, n as usize, sc.ge) {
            rr -= 1;
            trace.push(b'v');
        } else if cc > 0 && v == dp.get(rr, cc - 1) + horizontal_gap(ac, rr as usize, m as usize, sc.ge) {
            cc -= 1;
            trace.push(b'h');
        } else {
            if rr == 0 || cc == 0 {
                return false;
            }
            rr -= 1;
            cc -= 1;
            trace.push(b's');
        }
    }
    true
}

/// longNeedle: split alignment of s1 (consensus) against s2 (reference)
pub fn long_needle(s1: &[u8], s2: &[u8], ac: AlignConfig, sc: &DnaScore) -> Option<Align> {
    let m = s1.len() as i32;
    let n = s2.len() as i32;
    if m == 0 || n == 0 {
        return None;
    }
    let mut srev1 = s1.to_vec();
    reverse_complement(&mut srev1);
    let mut srev2 = s2.to_vec();
    reverse_complement(&mut srev2);

    let gain_max = sc.match_.max(sc.mismatch).max(sc.ge.max(0));
    let pen_max = (-sc.match_).max(-sc.mismatch).max((-sc.ge).max(0));
    let mut can_prune = sc.ge <= 0;
    if can_prune && !s1.iter().all(|&c| is_acgtn(c)) {
        can_prune = false;
    }
    if can_prune && !s2.iter().all(|&c| is_acgtn(c)) {
        can_prune = false;
    }
    let full_slack = m * gain_max + (m + n) * pen_max;

    let mut fmat = SplitDPMatrix::default();
    let mut rmat = SplitDPMatrix::default();
    let mut pm: Vec<i32> = Vec::new();
    let mut slack: i32 = 16;
    loop {
        let prune = can_prune && slack < full_slack;
        let lower_bound = m * gain_max - slack;
        split_fill(s1, s2, ac, sc, lower_bound, gain_max, prune, &mut fmat);
        split_fill(&srev1, &srev2, ac, sc, lower_bound, gain_max, prune, &mut rmat);
        if !prune && fmat.get(m, n) != rmat.get(m, n) {
            return None;
        }
        let mut best_score = fmat.get(m, n);
        let mut cons_left: i32 = 0;
        let mut ref_left: i32 = 0;
        for row in 0..=m {
            if fmat.empty(row) || rmat.empty(m - row) {
                continue;
            }
            let rlo = rmat.lo[(m - row) as usize];
            let rhi = rmat.hi[(m - row) as usize];
            let roff = rmat.off[(m - row) as usize];
            pm.resize((rhi - rlo + 1) as usize, 0);
            pm[0] = rmat.val[roff];
            for k in 1..=(rhi - rlo) as usize {
                pm[k] = pm[k - 1].max(rmat.val[roff + k]);
            }
            let foff = fmat.off[row as usize];
            let flo = fmat.lo[row as usize];
            for col in flo..=fmat.hi[row as usize] {
                let mut right = n - col;
                if right < rlo {
                    break;
                }
                if right > rhi {
                    right = rhi;
                }
                let score = fmat.val[foff + (col - flo) as usize] + pm[(right - rlo) as usize];
                if score > best_score {
                    best_score = score;
                    cons_left = row;
                    ref_left = col;
                }
            }
        }
        if prune && best_score < lower_bound {
            slack *= 4;
            if best_score > NEG / 2 && m * gain_max - best_score < slack {
                slack = m * gain_max - best_score;
            }
            continue;
        }
        if best_score == fmat.get(m, n) {
            return None;
        }
        let cons_right = m - cons_left;
        let mut ref_right: i32 = 0;
        let base = fmat.get(cons_left, ref_left);
        let r_end = rmat.hi[cons_right as usize].min(n - ref_left);
        let mut right = rmat.lo[cons_right as usize];
        while right <= r_end {
            if base + rmat.get(cons_right, right) == best_score {
                ref_right = right;
            }
            right += 1;
        }
        let mut trace: Vec<u8> = Vec::new();
        let mut rtrace: Vec<u8> = Vec::new();
        if !split_traceback(&fmat, cons_left, ref_left, m, n, ac, sc, &mut trace) || !split_traceback(&rmat, cons_right, ref_right, m, n, ac, sc, &mut rtrace) {
            if !prune {
                return None;
            }
            slack *= 4;
            continue;
        }
        let fwd = create_alignment_str(&trace, &s1[..cons_left as usize], &s2[..ref_left as usize]);
        let rvs = create_alignment_str(&rtrace, &srev1[..cons_right as usize], &srev2[..ref_right as usize]);
        let gapref = ((n - ref_right) - ref_left) as usize;
        let alilen = fwd.cols + rvs.cols + gapref;
        let mut al = Align::new(2, alilen);
        let j_end = rvs.cols;
        for i in 0..fwd.rows {
            let mut alicol = 0usize;
            while alicol < fwd.cols {
                al.set(i, alicol, fwd.get(i, alicol));
                alicol += 1;
            }
            for j in ref_left..(n - ref_right) {
                if i == 0 {
                    al.set(i, alicol, b'-');
                } else {
                    al.set(i, alicol, s2[j as usize]);
                }
                alicol += 1;
            }
            for j in 0..j_end {
                match rvs.get(i, j_end - j - 1) {
                    b'A' => al.set(i, alicol, b'T'),
                    b'C' => al.set(i, alicol, b'G'),
                    b'G' => al.set(i, alicol, b'C'),
                    b'T' => al.set(i, alicol, b'A'),
                    b'N' => al.set(i, alicol, b'N'),
                    b'-' => al.set(i, alicol, b'-'),
                    _ => {}
                }
                alicol += 1;
            }
        }
        return Some(al);
    }
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;

    #[test]
    fn hardware_fma_retains_fallback_alignment_bytes() {
        if !std::is_x86_feature_detected!("fma") { return; }
        let mut seed = 1u64;
        for case in 0..24 {
            let mut a = Align::new(1 + case % 7, 30 + case);
            let mut b = Align::new(1 + case % 5, 40 + case);
            for ch in a.data.iter_mut().chain(b.data.iter_mut()) {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                *ch = b"ACGTN-"[((seed >> 32) % 6) as usize];
            }
            for horizontal in [false, true] {
                for vertical in [false, true] {
                    let ac = AlignConfig { horizontal, vertical };
                    let sc = DnaScore::default();
                    let fallback = gotoh_impl(&a, &b, ac, &sc);
                    let accelerated = unsafe { gotoh_fma(&a, &b, ac, &sc) };
                    assert_eq!((accelerated.rows, accelerated.cols), (fallback.rows, fallback.cols));
                    assert_eq!(accelerated.data, fallback.data, "case {}", case);
                }
            }
        }
    }
}
