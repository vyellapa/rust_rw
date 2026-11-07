//! Port of cluster.h: split-read and paired-end clique clustering.

use crate::config::Config;
use crate::tags::*;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Edge {
    // ordering: weight, source, target (EdgeRecord::operator< is a total order)
    weight: u32,
    source: u32,
    target: u32,
}

type CompEdges = BTreeMap<u32, Vec<Edge>>;

// ---------------------------------------------------------------------------
// Split-read clustering
// ---------------------------------------------------------------------------

fn search_cliques_sr(c: &Config, comp_edge: &mut CompEdges, br: &mut [SRBamRecord], sv: &mut Vec<SVRecord>, svt: i32) {
    for edges in comp_edge.values_mut() {
        edges.sort_unstable();
        let first = edges[0];
        let mut clique: BTreeSet<u32> = BTreeSet::new();
        let mut incompatible: BTreeSet<u32> = BTreeSet::new();
        let mut seeds: BTreeSet<u64> = BTreeSet::new();
        let s = first.source as usize;
        clique.insert(first.source);
        seeds.insert(br[s].id);
        let chr = br[s].chr;
        let chr2 = br[s].chr2;
        let mut ciposlow = br[s].pos;
        let mut pos: u64 = br[s].pos as i64 as u64;
        let mut ciposhigh = br[s].pos;
        let mut ciendlow = br[s].pos2;
        let mut pos2: u64 = br[s].pos2 as i64 as u64;
        let mut ciendhigh = br[s].pos2;
        let mut mapq: i32 = br[s].qual;
        let mut inslen: i32 = br[s].inslen;

        let mut wiggle: u32 = c.max_read_sep;
        if is_translocation(svt) {
            wiggle = 2 * c.max_read_sep;
        } else {
            let mut svvar: u32 = (0.1 * (pos2.wrapping_sub(pos) as f64)).abs() as u32;
            if svt == 4 {
                svvar = (0.1 * inslen as f64).abs() as u32;
            }
            if wiggle < svvar {
                wiggle = svvar;
            }
            if wiggle > 1000 {
                wiggle = 1000;
            }
        }

        let mut grow = true;
        while grow {
            grow = false;
            for e in edges.iter() {
                if grow {
                    break;
                }
                let v: u32 = if !clique.contains(&e.source) && clique.contains(&e.target) {
                    e.source
                } else if clique.contains(&e.source) && !clique.contains(&e.target) {
                    e.target
                } else {
                    continue;
                };
                if incompatible.contains(&v) {
                    continue;
                }
                let bv = br[v as usize];
                if seeds.contains(&bv.id) {
                    continue;
                }
                let n_cpl = bv.pos.min(ciposlow);
                let n_cph = bv.pos.max(ciposhigh);
                let n_cel = bv.pos2.min(ciendlow);
                let n_ceh = bv.pos2.max(ciendhigh);
                if (n_cph - n_cpl) < wiggle as i32 && (n_ceh - n_cel) < wiggle as i32 && (!is_translocation(svt) || bv.chr2 == chr2) {
                    grow = true;
                }
                if grow {
                    clique.insert(v);
                    seeds.insert(bv.id);
                    ciposlow = n_cpl;
                    pos = pos.wrapping_add(bv.pos as i64 as u64);
                    ciposhigh = n_cph;
                    ciendlow = n_cel;
                    pos2 = pos2.wrapping_add(bv.pos2 as i64 as u64);
                    ciendhigh = n_ceh;
                    mapq += bv.qual;
                    inslen += bv.inslen;
                } else {
                    incompatible.insert(v);
                }
            }
        }

        // Count paired-end fragments only once
        let mut prev_seed: u64 = 0;
        let mut cliq_size: u32 = 0;
        let mut selected: BTreeSet<u64> = BTreeSet::new();
        for &sd in &seeds {
            if prev_seed.wrapping_add(1) != sd {
                selected.insert(sd);
                cliq_size += 1;
            }
            prev_seed = sd;
        }

        if cliq_size >= c.min_clique_size {
            let csz = clique.len() as u64;
            let sv_start = (pos / csz) as i32;
            let sv_end = (pos2 / csz) as i32;
            let sv_ins_len = inslen / clique.len() as i32;
            if sv_size_check_ins(sv_start, sv_end, svt, sv_ins_len) {
                if ciposlow > sv_start || ciposhigh < sv_start || ciendlow > sv_end || ciendhigh < sv_end {
                    eprintln!("Warning: Confidence intervals out of bounds: {},{},{}:{},{},{}", ciposlow, sv_start, ciposhigh, ciendlow, sv_end, ciendhigh);
                }
                let svid = sv.len() as i32;
                sv.push(SVRecord::new_sr(
                    chr,
                    sv_start,
                    chr2,
                    sv_end,
                    ciposlow - sv_start,
                    ciposhigh - sv_start,
                    ciendlow - sv_end,
                    ciendhigh - sv_end,
                    cliq_size as i32,
                    ((mapq as i64 as u64) / csz) as i32,
                    mapq,
                    sv_ins_len,
                    svt,
                    svid,
                ));
                for &m in &clique {
                    if selected.contains(&br[m as usize].id) {
                        br[m as usize].svid = svid;
                    }
                }
            }
        }
    }
}

/// cluster(c, vector<SRBamRecord>&, svs, svt); br must be sorted
pub fn cluster_sr(c: &Config, br: &mut [SRBamRecord], sv: &mut Vec<SVRecord>, svt: i32) {
    let mut comp: Vec<u32> = vec![0; br.len()];
    for ref_idx in 0..c.nchr {
        let br_lo = br.partition_point(|r| r.chr < ref_idx);
        let br_hi = br.partition_point(|r| r.chr <= ref_idx);
        if br_lo >= br_hi {
            continue;
        }
        let mut num_comp: u32 = 0;
        let mut comp_edge: CompEdges = BTreeMap::new();
        let mut last_connected = br_lo;
        let mut last_connected_start = br_lo;
        for i in br_lo..br_hi {
            if br[i].chr != ref_idx {
                continue;
            }
            if i > last_connected && !comp_edge.is_empty() {
                search_cliques_sr(c, &mut comp_edge, br, sv, svt);
                last_connected_start = last_connected;
                comp_edge.clear();
            }
            let mut varisize: u32 = c.max_read_sep;
            if is_translocation(svt) {
                varisize = 2 * c.max_read_sep;
            } else {
                let mut svvar: u32 = (0.1 * (br[i].pos2 - br[i].pos) as f64).abs() as u32;
                if svt == 4 {
                    svvar = (0.1 * br[i].inslen as f64).abs() as u32;
                }
                if varisize < svvar {
                    varisize = svvar;
                }
                if varisize > 1000 {
                    varisize = 1000;
                }
            }
            for j in (i + 1)..br_hi {
                if br[j].chr != ref_idx {
                    continue;
                }
                if (br[j].pos.wrapping_sub(br[i].pos)) as u32 > varisize {
                    break;
                }
                if svt == 4 && (br[j].inslen - br[i].inslen).unsigned_abs() > varisize {
                    continue;
                }
                if is_translocation(svt) && br[j].chr2 != br[i].chr2 {
                    continue;
                }
                if (br[j].pos2 - br[i].pos2).unsigned_abs() < varisize {
                    if j > last_connected {
                        last_connected = j;
                    }
                    let comp_index: u32;
                    if comp[i] == 0 {
                        if comp[j] == 0 {
                            num_comp += 1;
                            comp_index = num_comp;
                            comp[i] = comp_index;
                            comp[j] = comp_index;
                            comp_edge.entry(comp_index).or_default();
                        } else {
                            comp_index = comp[j];
                            comp[i] = comp_index;
                        }
                    } else if comp[j] == 0 {
                        comp_index = comp[i];
                        comp[j] = comp_index;
                    } else if comp[j] == comp[i] {
                        comp_index = comp[j];
                    } else {
                        let (ci, other) = if comp[j] < comp[i] { (comp[j], comp[i]) } else { (comp[i], comp[j]) };
                        comp_index = ci;
                        for k in last_connected_start..=last_connected {
                            if comp[k] == other {
                                comp[k] = ci;
                            }
                        }
                        let moved = comp_edge.remove(&other).unwrap_or_default();
                        comp_edge.get_mut(&ci).unwrap().extend(moved);
                    }
                    let list = comp_edge.get_mut(&comp_index).unwrap();
                    if (list.len() as u32) < c.graph_pruning {
                        let weight = (br[j].pos2 - br[i].pos2).unsigned_abs().wrapping_add((br[j].pos - br[i].pos).unsigned_abs());
                        list.push(Edge { weight, source: i as u32, target: j as u32 });
                    }
                }
            }
        }
        if !comp_edge.is_empty() {
            search_cliques_sr(c, &mut comp_edge, br, sv, svt);
            comp_edge.clear();
        }
        for k in br_lo..br_hi {
            comp[k] = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// Paired-end clustering
// ---------------------------------------------------------------------------

fn init_clique(el: &BamAlignRecord, svt: i32) -> (i32, i32, i32) {
    let (mut s, mut e, mut w) = (-1, -1, 0);
    if is_translocation(svt) {
        let ct = span_orientation(svt);
        if ct % 2 == 0 {
            s = el.pos + el.alen;
            e = if ct >= 2 { el.mpos } else { el.mpos + el.malen };
        } else {
            s = el.pos;
            e = if ct >= 2 { el.mpos + el.malen } else { el.mpos };
        }
        w = el.max_normal_isize;
    } else if svt == 0 {
        s = el.mpos + el.malen;
        e = el.pos + el.alen;
        w = el.max_normal_isize - el.alen.max(el.malen);
    } else if svt == 1 {
        s = el.mpos;
        e = el.pos;
        w = el.max_normal_isize - el.alen.max(el.malen);
    } else if svt == 2 {
        s = el.mpos + el.malen;
        e = el.pos;
        w = -el.max_normal_isize;
    } else if svt == 3 {
        s = el.mpos;
        e = el.pos + el.alen;
        w = el.max_normal_isize;
    }
    (s, e, w)
}

fn update_clique(el: &BamAlignRecord, sv_start: &mut i32, sv_end: &mut i32, wiggle: &mut i32, svt: i32) -> bool {
    if is_translocation(svt) {
        let ct = span_orientation(svt);
        let nss;
        let nse;
        let mut nw = *wiggle;
        if ct % 2 == 0 {
            nss = (*sv_start).max(el.pos + el.alen);
            nw -= nss - *sv_start;
            if ct >= 2 {
                nse = (*sv_end).min(el.mpos);
                nw -= *sv_end - nse;
            } else {
                nse = (*sv_end).max(el.mpos + el.malen);
                nw -= nse - *sv_end;
            }
        } else {
            nss = (*sv_start).min(el.pos);
            nw -= *sv_start - nss;
            if ct >= 2 {
                nse = (*sv_end).max(el.mpos + el.malen);
                nw -= nse - *sv_end;
            } else {
                nse = (*sv_end).min(el.mpos);
                nw -= *sv_end - nse;
            }
        }
        if nw > 0 {
            *sv_start = nss;
            *sv_end = nse;
            *wiggle = nw;
            return true;
        }
        false
    } else if svt == 0 || svt == 1 {
        let ct = span_orientation(svt);
        let nss;
        let nse;
        let mut nw;
        let wc;
        if ct == 0 {
            nss = (*sv_start).max(el.mpos + el.malen);
            nse = (*sv_end).max(el.pos + el.alen);
            nw = (el.max_normal_isize - (nss - el.mpos)).min(el.max_normal_isize - (nse - el.pos));
            wc = *wiggle - (nss - *sv_start).max(nse - *sv_end);
        } else {
            nss = (*sv_start).min(el.mpos);
            nse = (*sv_end).min(el.pos);
            nw = (el.max_normal_isize - (el.mpos + el.malen - nss)).min(el.max_normal_isize - (el.pos + el.alen - nse));
            wc = *wiggle - (*sv_start - nss).max(*sv_end - nse);
        }
        if wc < nw {
            nw = wc;
        }
        if nss < nse && nw >= 0 {
            *sv_start = nss;
            *sv_end = nse;
            *wiggle = nw;
            return true;
        }
        false
    } else if svt == 2 {
        let nss = (*sv_start).max(el.mpos + el.malen);
        let nse = (*sv_end).min(el.pos);
        let mut nw = el.pos + el.alen - el.mpos - el.max_normal_isize - (nse - nss);
        let wc = *wiggle + (*sv_end - *sv_start) - (nse - nss);
        if wc > nw {
            nw = wc;
        }
        if nss < nse && nw <= 0 {
            *sv_start = nss;
            *sv_end = nse;
            *wiggle = nw;
            return true;
        }
        false
    } else if svt == 3 {
        let nss = (*sv_start).min(el.mpos);
        let nse = (*sv_end).max(el.pos + el.alen);
        let mut nw = el.pos - (el.mpos + el.malen) + el.max_normal_isize - (nse - nss);
        let wc = *wiggle - ((nse - nss) - (*sv_end - *sv_start));
        if wc < nw {
            nw = wc;
        }
        if nss < nse && nw >= 0 {
            *sv_start = nss;
            *sv_end = nse;
            *wiggle = nw;
            return true;
        }
        false
    } else {
        false
    }
}

fn search_cliques_pe(c: &Config, comp_edge: &mut CompEdges, rec: &[BamAlignRecord], svs: &mut Vec<SVRecord>, svt: i32) {
    for edges in comp_edge.values_mut() {
        edges.sort_unstable();
        let first = edges[0];
        let s0 = &rec[first.source as usize];
        let mut clique: BTreeSet<u32> = BTreeSet::new();
        let mut incompatible: BTreeSet<u32> = BTreeSet::new();
        let cluster_ref = s0.tid;
        let cluster_mate_ref = s0.mtid;
        let (mut sv_start, mut sv_end, mut wiggle) = init_clique(s0, svt);
        if cluster_ref == cluster_mate_ref && sv_start >= sv_end {
            continue;
        }
        clique.insert(first.source);
        let mut grow = true;
        while grow {
            grow = false;
            for e in edges.iter() {
                if grow {
                    break;
                }
                let v = if !clique.contains(&e.source) && clique.contains(&e.target) {
                    e.source
                } else if clique.contains(&e.source) && !clique.contains(&e.target) {
                    e.target
                } else {
                    continue;
                };
                if incompatible.contains(&v) {
                    continue;
                }
                grow = update_clique(&rec[v as usize], &mut sv_start, &mut sv_end, &mut wiggle, svt);
                if grow {
                    clique.insert(v);
                } else {
                    incompatible.insert(v);
                }
            }
        }
        if clique.len() as u32 >= c.min_clique_size && sv_size_check(sv_start, sv_end, svt) {
            let mut r = SVRecord { chr: cluster_ref, chr2: cluster_mate_ref, ..Default::default() };
            if svt < DELLY_SVT_TRANS {
                r.sv_start = sv_start;
                r.sv_end = sv_end;
            } else {
                r.sv_start = (sv_start as u32).wrapping_add(1) as i32;
                r.sv_end = (sv_end as u32).wrapping_add(1) as i32;
            }
            r.pe_support = clique.len() as i32;
            let ci = wiggle.abs().max(50);
            r.ciposlow = -ci;
            r.ciposhigh = ci;
            r.ciendlow = -ci;
            r.ciendhigh = ci;
            r.mapq = 0;
            let mut mapqv: Vec<u8> = Vec::with_capacity(clique.len());
            for &m in &clique {
                mapqv.push(rec[m as usize].map_quality);
                r.mapq += rec[m as usize].map_quality as i32;
            }
            mapqv.sort_unstable();
            r.pe_map_quality = mapqv[mapqv.len() / 2] as i32;
            r.sr_support = 0;
            r.sr_align_quality = 0.0;
            r.precise = false;
            r.svt = svt;
            r.ins_len = 0;
            r.hom_len = 0;
            svs.push(r);
        }
    }
}

/// cluster(c, vector<BamAlignRecord>&, svs, varisize, svt); records must be sorted
pub fn cluster_pe(c: &Config, rec: &[BamAlignRecord], svs: &mut Vec<SVRecord>, varisize: u32, svt: i32) {
    let mut comp: Vec<u32> = vec![0; rec.len()];
    let mut num_comp: u32 = 0;
    let mut comp_edge: CompEdges = BTreeMap::new();
    let mut last_connected: usize = 0;
    let mut last_connected_start: usize = 0;
    for idx in 0..rec.len() {
        if idx > last_connected && !comp_edge.is_empty() {
            search_cliques_pe(c, &mut comp_edge, rec, svs, svt);
            last_connected_start = last_connected;
            comp_edge.clear();
        }
        let a = &rec[idx];
        let min_c = min_coord(a.pos, a.mpos, svt);
        let max_c = max_coord(a.pos, a.mpos, svt);
        let mut nidx = idx + 1;
        while nidx < rec.len() {
            let b = &rec[nidx];
            if (min_coord(b.pos, b.mpos, svt) + b.alen - min_c).unsigned_abs() > varisize {
                break;
            }
            if a.mtid != b.mtid {
                nidx += 1;
                continue;
            }
            let bmin = min_coord(b.pos, b.mpos, svt);
            let bmax = max_coord(b.pos, b.mpos, svt);
            if pairs_disagree(min_c, max_c, a.alen, a.max_normal_isize, bmin, bmax, b.alen, b.max_normal_isize, svt) {
                nidx += 1;
                continue;
            }
            if nidx > last_connected {
                last_connected = nidx;
            }
            let comp_index: u32;
            if comp[idx] == 0 {
                if comp[nidx] == 0 {
                    num_comp += 1;
                    comp_index = num_comp;
                    comp[idx] = comp_index;
                    comp[nidx] = comp_index;
                    comp_edge.entry(comp_index).or_default();
                } else {
                    comp_index = comp[nidx];
                    comp[idx] = comp_index;
                }
            } else if comp[nidx] == 0 {
                comp_index = comp[idx];
                comp[nidx] = comp_index;
            } else if comp[nidx] == comp[idx] {
                comp_index = comp[nidx];
            } else {
                let (ci, other) = if comp[nidx] < comp[idx] { (comp[nidx], comp[idx]) } else { (comp[idx], comp[nidx]) };
                comp_index = ci;
                for k in last_connected_start..=last_connected {
                    if comp[k] == other {
                        comp[k] = ci;
                    }
                }
                let moved = comp_edge.remove(&other).unwrap_or_default();
                comp_edge.get_mut(&ci).unwrap().extend(moved);
            }
            let list = comp_edge.get_mut(&comp_index).unwrap();
            if (list.len() as u32) < c.graph_pruning {
                let d = ((bmin - min_c) - (bmax - max_c)).abs() - (a.median - b.median).abs();
                let weight = ((d.abs() as f64 + 1.0).log2()) as u32;
                list.push(Edge { weight, source: idx as u32, target: nidx as u32 });
            }
            nidx += 1;
        }
    }
    if !comp_edge.is_empty() {
        search_cliques_pe(c, &mut comp_edge, rec, svs, svt);
        comp_edge.clear();
    }
}
