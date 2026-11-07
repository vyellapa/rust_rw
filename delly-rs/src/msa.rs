//! Port of msa.h: UPGMA guide tree + progressive profile alignment + consensus.

use crate::align::{gotoh, Align, AlignConfig};
use crate::config::Config;

/// Longest common subsequence length
pub fn lcs(s1: &[u8], s2: &[u8]) -> i32 {
    let m = s1.len();
    let n = s2.len();
    let mut prevdiag: i32;
    let mut onecol = vec![0i32; n + 1];
    // row 0 is all zeros
    for i in 1..=m {
        prevdiag = 0;
        onecol[0] = 0;
        let c1 = s1[i - 1];
        for j in 1..=n {
            let prevprevdiag = prevdiag;
            prevdiag = onecol[j];
            if c1 == s2[j - 1] {
                onecol[j] = prevprevdiag + 1;
            } else if onecol[j - 1] > onecol[j] {
                onecol[j] = onecol[j - 1];
            }
        }
    }
    onecol[n]
}

struct Mat {
    n: usize,
    d: Vec<i64>,
}

impl Mat {
    #[inline]
    fn g(&self, i: usize, j: usize) -> i64 {
        self.d[i * self.n + j]
    }
    #[inline]
    fn s(&mut self, i: usize, j: usize, v: i64) {
        self.d[i * self.n + j] = v;
    }
}

fn closest_pair(d: &Mat, num: usize, di: &mut usize, dj: &mut usize) -> i64 {
    let mut dmax: i64 = -1;
    for i in 0..num {
        for j in (i + 1)..num {
            let v = d.g(i, j);
            if v > dmax {
                dmax = v;
                *di = i;
                *dj = j;
            }
        }
    }
    dmax
}

fn update_distance_matrix(d: &mut Mat, p: &[[i64; 3]], num: usize, di: usize, dj: usize) {
    for i in 0..num {
        if p[i][0] == -1 {
            let a = if di < i { d.g(di, i) } else { d.g(i, di) };
            let b = if dj < i { d.g(dj, i) } else { d.g(i, dj) };
            d.s(i, num, (a + b) / 2);
        }
    }
    for i in 0..di {
        d.s(i, di, -1);
    }
    for i in (di + 1)..(num + 1) {
        d.s(di, i, -1);
    }
    for i in 0..dj {
        d.s(i, dj, -1);
    }
    for i in (dj + 1)..(num + 1) {
        d.s(dj, i, -1);
    }
}

fn palign(c: &Config, sps: &[&[u8]], p: &[[i64; 3]], root: usize) -> Align {
    if p[root][1] == -1 && p[root][2] == -1 {
        Align::from_seq(sps[root])
    } else {
        let a1 = palign(c, sps, p, p[root][1] as usize);
        let a2 = palign(c, sps, p, p[root][2] as usize);
        gotoh(&a1, &a2, AlignConfig { horizontal: true, vertical: true }, &c.aliscore)
    }
}

/// consensus(c, align, cs)
fn consensus(c: &Config, al: &Align) -> Vec<u8> {
    let rows = al.rows;
    let cols = al.cols;
    let mut fl = vec![false; rows * cols];
    let mut cov = vec![0i32; cols];
    for i in 0..rows {
        let mut start: i64 = 0;
        let mut end: i64 = -1;
        for j in 0..cols {
            if al.get(i, j) != b'-' {
                end = j as i64;
            } else if end == -1 {
                start = j as i64 + 1;
            }
        }
        let mut j = start;
        while j <= end {
            cov[j as usize] += 1;
            fl[i * cols + j as usize] = true;
            j += 1;
        }
    }
    let cov_threshold = 2.max((c.min_clique_size as i32).min(rows as i32));
    let mut cons = vec![b'-'; cols];
    for j in 0..cols {
        let mut max_idx = 4;
        if cov[j] >= cov_threshold {
            let mut count = [0i32; 5];
            for i in 0..rows {
                if fl[i * cols + j] {
                    match al.get(i, j) {
                        b'A' | b'a' => count[0] += 1,
                        b'C' | b'c' => count[1] += 1,
                        b'G' | b'g' => count[2] += 1,
                        b'T' | b't' => count[3] += 1,
                        _ => count[4] += 1,
                    }
                }
            }
            max_idx = 0;
            let mut max_count = count[0];
            for (k, &v) in count.iter().enumerate().skip(1) {
                if v > max_count {
                    max_count = v;
                    max_idx = k;
                }
            }
        }
        match max_idx {
            0 => cons[j] = b'A',
            1 => cons[j] = b'C',
            2 => cons[j] = b'G',
            3 => cons[j] = b'T',
            _ => {}
        }
    }
    cons.into_iter().filter(|&ch| ch != b'-').collect()
}

/// msa: the input sequences must be in the iteration order of the C++ container.
pub fn msa(c: &Config, sps: &[&[u8]]) -> Vec<u8> {
    let num = sps.len();
    let dim = 2 * num + 1;
    let mut d = Mat { n: dim, d: vec![0i64; dim * dim] };
    for i in 0..dim {
        for j in (i + 1)..dim {
            d.s(i, j, -1);
        }
    }
    for i in 0..num {
        for j in (i + 1)..num {
            let l = lcs(sps[i], sps[j]) as i64;
            let mn = sps[i].len().min(sps[j].len()) as i64;
            // (int * 100) / size_t: performed in unsigned 64-bit arithmetic
            d.s(i, j, (((l * 100) as u64) / (mn as u64)) as i32 as i64);
        }
    }
    let mut p = vec![[-1i64; 3]; dim];
    let mut nn = num;
    while nn < 2 * num + 1 {
        let mut di = 0usize;
        let mut dj = 0usize;
        if closest_pair(&d, nn, &mut di, &mut dj) == -1 {
            break;
        }
        p[di][0] = nn as i64;
        p[dj][0] = nn as i64;
        p[nn][1] = di as i64;
        p[nn][2] = dj as i64;
        update_distance_matrix(&mut d, &p, nn, di, dj);
        nn += 1;
    }
    let root = if nn > 0 { nn - 1 } else { 0 };
    let al = palign(c, sps, &p, root);
    consensus(c, &al)
}
