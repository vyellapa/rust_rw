//! Port of bolog.h (genotype likelihoods) and modvcf.h (vcfParse / vcfOutput).

use crate::config::Config;
use crate::coverage::Annotation;
use crate::hts::*;
use crate::ploidy::sv_ploidy;
use crate::tags::*;
use crate::util::*;

const SMALLEST_GL: f64 = -1000.0;
/// clang constant-folds log10(2.0) and negates it; this is the exact double it emits.
const NEG_LOG10_2: u64 = 0xBFD34413509F79FF;
const DELLY_VERSION: &str = "2.8.1";

pub struct BoLog {
    phred2prob: Vec<f64>,
}

impl BoLog {
    pub fn new() -> Self {
        let n = (-10.0 * SMALLEST_GL).round() as i32;
        BoLog { phred2prob: (0..=n).map(|i| 10f64.powf(-(i as f64 / 10.0))).collect() }
    }
}

/// _computeGLs
#[allow(clippy::too_many_arguments)]
fn compute_gls(bl: &BoLog, mref: &[u8], malt: &[u8], gls: &mut [f32], gq: &mut [i32], gts: &mut [i32], f: usize, ploidy: u8) {
    let p = &bl.phred2prob;
    let mut gl = [0f64; 3];
    let mut pe_depth: u32 = (mref.len() + malt.len()) as u32;
    if ploidy == 0 {
        pe_depth = 0;
    }
    for &q in mref {
        let x = p[q as usize];
        gl[0] += x.log10();
        gl[1] += (x + (1.0 - x)).log10();
        gl[2] += (1.0 - x).log10();
    }
    for &q in malt {
        let x = p[q as usize];
        gl[0] += (1.0 - x).log10();
        gl[1] += ((1.0 - x) + x).log10();
        gl[2] += x.log10();
    }
    gl[1] = (pe_depth as f64).mul_add(f64::from_bits(NEG_LOG10_2), gl[1]);
    if ploidy == 1 {
        gl[1] = SMALLEST_GL + gl[0] + gl[2];
    }
    let mut best = 0usize;
    let mut best_val = gl[0];
    for g in 1..=2 {
        if gl[g] >= best_val {
            best_val = gl[g];
            best = g;
        }
    }
    for g in gl.iter_mut() {
        *g -= best_val;
        *g = if *g > SMALLEST_GL { *g } else { SMALLEST_GL };
    }
    let pl: [u32; 3] = [(-10.0 * gl[0]).round() as u32, (-10.0 * gl[1]).round() as u32, (-10.0 * gl[2]).round() as u32];
    if pe_depth != 0 && pl[0].wrapping_add(pl[1]).wrapping_add(pl[2]) > 0 {
        let probsum = if ploidy == 1 { p[pl[0] as usize] + p[pl[2] as usize] } else { p[pl[0] as usize] + p[pl[1] as usize] + p[pl[2] as usize] };
        let mut likelihood = (1.0 - 1.0 / probsum).log10();
        likelihood = if likelihood > SMALLEST_GL { likelihood } else { SMALLEST_GL };
        gq[f] = (-10.0 * likelihood).round() as i32;
        if ploidy == 1 {
            gts[f * 2] = if best == 0 { bcf_gt_unphased(1) } else { bcf_gt_unphased(0) };
            gts[f * 2 + 1] = BCF_INT32_VECTOR_END;
        } else if best == 0 {
            gts[f * 2] = bcf_gt_unphased(1);
            gts[f * 2 + 1] = bcf_gt_unphased(1);
        } else if best == 1 {
            gts[f * 2] = bcf_gt_unphased(0);
            gts[f * 2 + 1] = bcf_gt_unphased(1);
        } else {
            gts[f * 2] = bcf_gt_unphased(0);
            gts[f * 2 + 1] = bcf_gt_unphased(0);
        }
    } else {
        gts[f * 2] = BCF_GT_MISSING;
        gts[f * 2 + 1] = if ploidy == 1 { BCF_INT32_VECTOR_END } else { BCF_GT_MISSING };
        gq[f] = 0;
    }
    gls[f * 3 + 2] = gl[0] as f32;
    gls[f * 3 + 1] = gl[1] as f32;
    gls[f * 3] = gl[2] as f32;
}

/// _locusGLs
fn locus_gls(jc: &crate::coverage::JunctionCount, gls: &mut [f32], gq: &mut [i32], gts: &mut [i32], f: usize, ploidy: u8) {
    gls[f * 3] = jc.gl[0];
    gls[f * 3 + 1] = jc.gl[1];
    gls[f * 3 + 2] = jc.gl[2];
    gq[f] = jc.locus_gq;
    if ploidy == 0 {
        gts[f * 2] = BCF_GT_MISSING;
        gts[f * 2 + 1] = BCF_GT_MISSING;
        gq[f] = 0;
    } else if ploidy == 1 {
        gts[f * 2] = if jc.locus_gt > 0 { bcf_gt_unphased(1) } else { bcf_gt_unphased(0) };
        gts[f * 2 + 1] = BCF_INT32_VECTOR_END;
    } else {
        gts[f * 2] = if jc.locus_gt == 2 { bcf_gt_unphased(1) } else { bcf_gt_unphased(0) };
        gts[f * 2 + 1] = if jc.locus_gt >= 1 { bcf_gt_unphased(1) } else { bcf_gt_unphased(0) };
    }
}

/// _replaceIUPAC
fn replace_iupac(alleles: &str) -> String {
    let a = alleles.as_bytes();
    let mut out = vec![0u8; a.len()];
    let mut in_tag = 0;
    let mut in_ref = true;
    for (i, &ch) in a.iter().enumerate() {
        if ch == b',' {
            in_ref = false;
        }
        if in_ref || in_tag != 0 || matches!(ch, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' | b'<' | b'>' | b']' | b'[' | b',') {
            out[i] = ch;
            if ch == b'<' {
                in_tag = 1;
            } else if ch == b']' {
                in_tag = 2;
            } else if ch == b'[' {
                in_tag = 3;
            } else if ch == b'>' && in_tag == 1 {
                in_tag = 0;
            } else if ch == b']' && in_tag == 2 {
                in_tag = 0;
            } else if ch == b'[' && in_tag == 3 {
                in_tag = 0;
            }
        } else {
            out[i] = match ch {
                b'U' | b'u' => b'T',
                b'R' | b'r' => b'A',
                b'Y' | b'y' => b'C',
                b'S' | b's' => b'C',
                b'W' | b'w' => b'A',
                b'K' | b'k' => b'G',
                b'M' | b'm' => b'A',
                b'B' | b'b' => b'C',
                b'D' | b'd' => b'A',
                b'H' | b'h' => b'A',
                b'V' | b'v' => b'A',
                _ => b'N',
            };
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn is_symbolic(alleles: &str) -> (bool, Option<usize>) {
    match alleles.find(',') {
        None => (true, None),
        Some(cp) => {
            let alt = &alleles[cp + 1..];
            let sym = !alt.is_empty() && (alt.starts_with('<') || alt.contains('[') || alt.contains(']'));
            (sym, Some(cp))
        }
    }
}

/// vcfOutput
pub fn vcf_output(c: &Config, hdr: &Header, svs: &[SVRecord], ann: &Annotation) {
    let bl = BoLog::new();
    let fmtout = if c.outfile == "-" { "w" } else { "wb" };
    let mut w = match BcfWriter::open(&c.outfile, fmtout) {
        Some(w) => w,
        None => {
            eprintln!("Fail to open file {}", c.outfile);
            return;
        }
    };
    w.hdr_append(&format!("##fileDate={}", today_iso()));
    for line in [
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##ALT=<ID=BND,Description=\"Translocation\">",
        "##ALT=<ID=INS,Description=\"Insertion\">",
        "##FILTER=<ID=LowQual,Description=\"Poor quality and insufficient number of PEs and SRs.\">",
        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">",
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">",
        "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 coordinate in case of an inter-chromosomal translocation\">",
        "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal translocation\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">",
        "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">",
        "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">",
        "##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">",
        "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">",
        "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">",
        "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">",
        "##INFO=<ID=CONSBP,Number=1,Type=Integer,Description=\"Consensus SV breakpoint position\">",
        "##INFO=<ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">",
        "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length; negative for DEL, positive for DUP/INV/INS.\">",
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",
        "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">",
        "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">",
        "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Breakpoint homology length\">",
        "##INFO=<ID=SUBTYPE,Number=1,Type=String,Description=\"SV subtype: INS:ME:ALU, INS:ME:LINE1, INS:ME:SVA, INS:NUMT, INS:LTR, INS:HERVK, INS:TR, or DEL:TR\">",
        "##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description=\"Identifier of the merged locus\">",
        "##INFO=<ID=NALLELE,Number=1,Type=Integer,Description=\"Number of distinct alleles at this locus\">",
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">",
        "##INFO=<ID=INSSTRAND,Number=1,Type=String,Description=\"Insertion strand for MEIs\">",
        "##INFO=<ID=TRPERIOD,Number=1,Type=Integer,Description=\"Tandem repeat period in bp\">",
        "##INFO=<ID=TRSTART,Number=1,Type=Integer,Description=\"Start of the reference tandem repeat\">",
        "##INFO=<ID=TREND,Number=1,Type=Integer,Description=\"End of the reference tandem repeat\">",
        "##INFO=<ID=TRCOPIES,Number=1,Type=Float,Description=\"Tandem repeat copy number\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">",
        "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the SV\">",
        "##FORMAT=<ID=RCL,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the left control region\">",
        "##FORMAT=<ID=RCR,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the right control region\">",
        "##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description=\"Read-depth based copy-number estimate for autosomal sites\">",
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">",
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">",
        "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">",
        "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">",
        "##FORMAT=<ID=HP,Number=4,Type=Integer,Description=\"Haplotype-specific junction read counts (HP1_ref,HP1_alt,HP2_ref,HP2_alt)\">",
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier from HP-tagged alignments\">",
        "##FORMAT=<ID=MR,Number=4,Type=Integer,Description=\"Methylation % for REF allele [SV start left/right, SV end left/right]\">",
        "##FORMAT=<ID=MA,Number=4,Type=Integer,Description=\"Methylation % for ALT allele [SV start left/right, SV end left/right]\">",
        "##FORMAT=<ID=MNC,Number=4,Type=Integer,Description=\"Unique CpG sites observed per window [SV start left/right, SV end left/right]\">",
        "##FORMAT=<ID=MDV,Number=4,Type=Integer,Description=\"Avg. read depth per CpG site per window [SV start left/right, SV end left/right]\">",
    ] {
        w.hdr_append(line);
    }
    w.hdr_append(&format!("##reference={}", c.genome));
    for i in 0..hdr.names.len() {
        w.hdr_append(&format!("##contig=<ID={},length={}>", hdr.names[i], hdr.lens[i]));
    }
    let nfiles = ann.jct.len();
    for f in 0..nfiles {
        w.add_sample(Some(&c.sample_name[f]));
    }
    w.add_sample(None);
    if !w.write_header() {
        eprintln!("Error: Failed to write BCF header!");
    }

    if !svs.is_empty() {
        let ns = w.nsamples() as usize;
        let mut gts = vec![0i32; ns * 2];
        let mut gls = vec![0f32; ns * 3];
        let mut rcl = vec![0i32; ns];
        let mut rc = vec![0i32; ns];
        let mut rcr = vec![0i32; ns];
        let mut cnest = vec![0i32; ns];
        let mut drc = vec![0i32; ns];
        let mut dvc = vec![0i32; ns];
        let mut rrc = vec![0i32; ns];
        let mut rvc = vec![0i32; ns];
        let mut gq = vec![0i32; ns];
        let mut hp = vec![0i32; ns * 4];
        let mut ps = vec![0i32; ns];
        let mut ma = vec![0i32; ns * 4];
        let mut mr = vec![0i32; ns * 4];
        let mut mnc = vec![0i32; ns * 4];
        let mut mdv = vec![0i32; ns * 4];
        let mut pl = vec![0i32; ns * 3];
        let mut ft: Vec<String> = vec![String::new(); ns];
        log("Genotyping");

        // Multi-allelic loci (only relevant for re-genotyped merged sites)
        let mut locus_members: std::collections::BTreeMap<i32, Vec<i32>> = std::collections::BTreeMap::new();
        let mut gt_class: std::collections::BTreeMap<i32, Vec<u8>> = std::collections::BTreeMap::new();
        for sv in svs {
            if sv.alleleid < 0 || sv.nallele < 2 {
                continue;
            }
            if sv.sr_support == 0 && sv.pe_support == 0 {
                continue;
            }
            let id = sv.id as usize;
            if !(0..nfiles).any(|f| ann.jct[f][id].joint) {
                continue;
            }
            locus_members.entry(sv.alleleid).or_default().push(sv.id);
            let mut cls = vec![3u8; nfiles];
            for f in 0..nfiles {
                let mut ploidy = 2u8;
                if f < c.sex_model.sex.len() {
                    ploidy = sv_ploidy(&c.sex_model, c.sex_model.sex[f], sv.chr, sv.sv_start, sv.chr2, sv.sv_end);
                }
                if ploidy != 2 || !ann.jct[f][id].joint {
                    continue;
                }
                if sv.precise && ann.jct[f][id].locus_gt >= 0 {
                    locus_gls(&ann.jct[f][id], &mut gls, &mut gq, &mut gts, f, ploidy);
                } else if sv.precise {
                    compute_gls(&bl, &ann.jct[f][id].ref_, &ann.jct[f][id].alt, &mut gls, &mut gq, &mut gts, f, ploidy);
                } else {
                    compute_gls(&bl, &ann.span[f][id].ref_, &ann.span[f][id].alt, &mut gls, &mut gq, &mut gts, f, ploidy);
                }
                if gts[f * 2] == BCF_GT_MISSING {
                    continue;
                }
                cls[f] = (bcf_gt_allele(gts[f * 2]) + bcf_gt_allele(gts[f * 2 + 1])) as u8;
            }
            gt_class.insert(sv.id, cls);
        }

        w.init_rec();
        for sv in svs {
            if sv.sr_support == 0 && sv.pe_support == 0 {
                continue;
            }
            let id = sv.id as usize;
            if !c.has_vcf_file {
                let mut tot: u32 = 0;
                for f in 0..nfiles {
                    tot += (ann.span[f][id].alt.len() + ann.jct[f][id].alt.len()) as u32;
                }
                if tot < 2 {
                    continue;
                }
            }
            let mut tmpi = w.id2int("PASS");
            if sv.chr == sv.chr2 {
                if (sv.pe_support < 3 || sv.pe_map_quality < 20) && (sv.sr_support < 3 || sv.sr_map_quality < 20) {
                    tmpi = w.id2int("LowQual");
                }
            } else if (sv.pe_support < 5 || sv.pe_map_quality < 20) && (sv.sr_support < 5 || sv.sr_map_quality < 20) {
                tmpi = w.id2int("LowQual");
            }
            let rid = w.name2id(&hdr.names[sv.chr as usize]);
            w.set_rid(rid);
            let sv_start_pos = (sv.sv_start - 1).max(0);
            let mut sv_end_pos = sv.sv_end;
            if sv_end_pos < 1 {
                sv_end_pos = 1;
            }
            if sv_end_pos > hdr.lens[sv.chr2 as usize] as i32 {
                sv_end_pos = hdr.lens[sv.chr2 as usize] as i32;
            }
            w.set_pos(sv_start_pos as i64);
            w.update_id(&format!("{}{:08}", add_id(sv.svt), sv.id));
            let alleles = replace_iupac(&sv.alleles);
            w.update_alleles_str(&alleles);
            w.update_filter(tmpi);
            if sv.precise {
                w.info_flag("PRECISE");
            } else {
                w.info_flag("IMPRECISE");
            }
            w.info_string("SVTYPE", add_id(sv.svt));
            w.info_string("SVMETHOD", &format!("EMBL.DELLYv{}", DELLY_VERSION));
            if sv.svt < DELLY_SVT_TRANS {
                let (sym, cp) = is_symbolic(&alleles);
                let v = if !sym {
                    sv_start_pos + cp.unwrap() as i32
                } else {
                    if sv_end_pos < sv_start_pos + 1 {
                        sv_end_pos = sv_start_pos + 1;
                    }
                    sv_end_pos
                };
                w.info_i32("END", &[v]);
            } else {
                w.info_i32("END", &[sv_start_pos + 1]);
                w.info_string("CHR2", &hdr.names[sv.chr2 as usize]);
                w.info_i32("POS2", &[sv_end_pos]);
            }
            if sv.svt == 4 {
                w.info_i32("SVLEN", &[sv.ins_len]);
            } else if sv.svt == 2 {
                w.info_i32("SVLEN", &[sv.sv_start - sv.sv_end]);
            } else if sv.svt == 3 || sv.svt == 0 || sv.svt == 1 {
                w.info_i32("SVLEN", &[sv.sv_end - sv.sv_start]);
            }
            w.info_i32("PE", &[sv.pe_support]);
            w.info_i32("MAPQ", &[sv.pe_map_quality]);
            w.info_string("CT", add_orientation(sv.svt));
            w.info_i32("CIPOS", &[sv.ciposlow, sv.ciposhigh]);
            w.info_i32("CIEND", &[sv.ciendlow, sv.ciendhigh]);
            if sv.alleleid >= 0 {
                w.info_i32("ALLELEID", &[sv.alleleid]);
                w.info_i32("NALLELE", &[sv.nallele]);
            }
            if sv.tr_start >= 0 {
                w.info_i32("TRSTART", &[sv.tr_start + 1]);
                w.info_i32("TREND", &[sv.tr_end]);
            }
            if sv.precise {
                w.info_i32("SRMAPQ", &[sv.sr_map_quality]);
                w.info_i32("INSLEN", &[sv.ins_len]);
                w.info_i32("HOMLEN", &[sv.hom_len]);
                w.info_i32("SR", &[sv.sr_support]);
                w.info_f32("SRQ", &[sv.sr_align_quality]);
                if !sv.consensus.is_empty() {
                    w.info_string("CONSENSUS", &sv.consensus);
                    w.info_f32("CE", &[entropy(sv.consensus.as_bytes()) as f32]);
                    w.info_i32("CONSBP", &[sv.cons_bp]);
                }
            }
            if !is_translocation(sv.svt) {
                if sv.anno.hom_len > 0 {
                    w.info_i32("HOMLEN", &[sv.anno.hom_len]);
                }
                if sv.anno.seq_type > 0 && sv.anno.seq_type < 7 {
                    const ST: [&str; 7] = ["", "INS:ME:ALU", "INS:ME:LINE1", "INS:ME:SVA", "INS:NUMT", "INS:LTR", "INS:HERVK"];
                    w.info_string("SUBTYPE", ST[sv.anno.seq_type as usize]);
                    w.info_string("INSSTRAND", if sv.anno.is_rc { "-" } else { "+" });
                } else if sv.anno.seq_type == 7 {
                    w.info_string("SUBTYPE", if sv.svt == 4 { "INS:TR" } else { "DEL:TR" });
                    w.info_i32("TRPERIOD", &[sv.anno.tr_period]);
                    w.info_f32("TRCOPIES", &[sv.anno.tr_copies]);
                } else if sv.tr_start >= 0 && sv.nallele > 1 && (sv.svt == 2 || sv.svt == 4) {
                    w.info_string("SUBTYPE", if sv.svt == 4 { "INS:TR" } else { "DEL:TR" });
                }
            }

            for f in 0..nfiles {
                let jc = &ann.jct[f][id];
                let sp = &ann.span[f][id];
                drc[f] = sp.ref_.len() as i32;
                dvc[f] = sp.alt.len() as i32;
                rrc[f] = jc.ref_.len() as i32;
                rvc[f] = jc.alt.len() as i32;
                hp[f * 4] = jc.hp1ref.len() as i32;
                hp[f * 4 + 1] = jc.hp1alt.len() as i32;
                hp[f * 4 + 2] = jc.hp2ref.len() as i32;
                hp[f * 4 + 3] = jc.hp2alt.len() as i32;
                ps[f] = jc.ps;
                for k in 0..4 {
                    ma[f * 4 + k] = BCF_INT32_MISSING;
                    mr[f * 4 + k] = BCF_INT32_MISSING;
                    mnc[f * 4 + k] = BCF_INT32_MISSING;
                    mdv[f * 4 + k] = BCF_INT32_MISSING;
                }
                let mut ploidy = 2u8;
                if f < c.sex_model.sex.len() {
                    ploidy = sv_ploidy(&c.sex_model, c.sex_model.sex[f], sv.chr, sv.sv_start, sv.chr2, sv.sv_end);
                }
                if sv.precise && jc.locus_gt >= 0 {
                    locus_gls(jc, &mut gls, &mut gq, &mut gts, f, ploidy);
                } else if sv.precise {
                    compute_gls(&bl, &jc.ref_, &jc.alt, &mut gls, &mut gq, &mut gts, f, ploidy);
                } else {
                    compute_gls(&bl, &sp.ref_, &sp.alt, &mut gls, &mut gq, &mut gts, f, ploidy);
                }
                if gts[f * 2] == BCF_GT_MISSING {
                    for k in 0..3 {
                        gls[f * 3 + k] = f32::from_bits(BCF_FLOAT_MISSING_BITS);
                        pl[f * 3 + k] = BCF_INT32_MISSING;
                    }
                } else if ploidy == 1 {
                    pl[f * 3] = 0.0f32.max((-10.0f32 * gls[f * 3]).round()) as i32;
                    pl[f * 3 + 1] = 0.0f32.max((-10.0f32 * gls[f * 3 + 2]).round()) as i32;
                    pl[f * 3 + 2] = BCF_INT32_VECTOR_END;
                } else {
                    for k in 0..3 {
                        pl[f * 3 + k] = 0.0f32.max((-10.0f32 * gls[f * 3 + k]).round()) as i32;
                    }
                }
                if ploidy == 2 && ps[f] != -1 {
                    let hp1a = hp[f * 4 + 1];
                    let hp2a = hp[f * 4 + 3];
                    let het = gts[f * 2] == bcf_gt_unphased(0) && gts[f * 2 + 1] == bcf_gt_unphased(1);
                    if het && (hp1a + hp2a) > 0 && hp1a != hp2a {
                        if hp1a > hp2a {
                            gts[f * 2] = bcf_gt_phased(1);
                            gts[f * 2 + 1] = bcf_gt_phased(0);
                        } else {
                            gts[f * 2] = bcf_gt_phased(0);
                            gts[f * 2 + 1] = bcf_gt_phased(1);
                        }
                    }
                }
                if ploidy == 2 && ps[f] == -1 && gt_class.contains_key(&sv.id) {
                    let het = gts[f * 2] == bcf_gt_unphased(0) && gts[f * 2 + 1] == bcf_gt_unphased(1);
                    if het {
                        let mem = locus_members.get(&sv.alleleid).cloned().unwrap_or_default();
                        let mut n_het = 0;
                        let mut n_other = 0;
                        let mut first_het = -1;
                        for &m in &mem {
                            let cl = gt_class.entry(m).or_default().get(f).copied().unwrap_or(0);
                            if cl == 1 {
                                n_het += 1;
                                if first_het < 0 {
                                    first_het = m;
                                }
                            } else if cl == 2 {
                                n_other += 1;
                            }
                        }
                        if n_het == 2 && n_other == 0 {
                            if sv.id == first_het {
                                gts[f * 2] = bcf_gt_phased(1);
                                gts[f * 2 + 1] = bcf_gt_phased(0);
                            } else {
                                gts[f * 2] = bcf_gt_phased(0);
                                gts[f * 2 + 1] = bcf_gt_phased(1);
                            }
                            ps[f] = sv.alleleid;
                        }
                    }
                }
                let r = &ann.rc[f][id];
                rcl[f] = r.left_rc;
                rc[f] = r.rc;
                rcr[f] = r.right_rc;
                cnest[f] = -1;
                if (rcl[f] + rcr[f]) > 0 && ploidy > 0 {
                    let mut cn = ploidy as f64 * rc[f] as f64 / (rcl[f] + rcr[f]) as f64;
                    if cn < 0.0 {
                        cn = 0.0;
                    }
                    if cn > 100000.0 {
                        cn = 100000.0;
                    }
                    cnest[f] = cn.round() as i32;
                }
                ft[f] = if gq[f] < 15 { "LowQual".to_string() } else { "PASS".to_string() };
            }
            let mut q = sv.mapq;
            if q < 0 {
                q = 0;
            }
            if q > 10000 {
                q = 10000;
            }
            w.set_qual(q as f32);
            let mut ac = 0;
            let mut an = 0;
            for &g in gts.iter().take(ns * 2) {
                if g == BCF_INT32_VECTOR_END || bcf_gt_is_missing(g) {
                    continue;
                }
                an += 1;
                if bcf_gt_allele(g) > 0 {
                    ac += 1;
                }
            }
            w.info_i32("AC", &[ac]);
            w.info_i32("AN", &[an]);
            w.genotypes(&gts);
            w.format_i32("GQ", &gq);
            w.format_i32("PL", &pl);
            w.format_string("FT", &ft);
            w.format_i32("RCL", &rcl);
            w.format_i32("RC", &rc);
            w.format_i32("RCR", &rcr);
            w.format_i32("RDCN", &cnest);
            w.format_i32("DR", &drc);
            w.format_i32("DV", &dvc);
            w.format_i32("RR", &rrc);
            w.format_i32("RV", &rvc);
            w.format_i32("HP", &hp);
            w.format_i32("PS", &ps);
            w.format_i32("MR", &mr);
            w.format_i32("MA", &ma);
            w.format_i32("MNC", &mnc);
            w.format_i32("MDV", &mdv);
            w.write_rec();
        }
    }
    w.close();
    if c.outfile != "-" {
        bcf_index_build_csi(&c.outfile, 14);
    }
}

/// vcfParse: load SV sites from a delly VCF/BCF for genotyping
pub fn vcf_parse(c: &Config, hd: &Header, svs: &mut Vec<SVRecord>) {
    let mut r = match BcfReader::open(&c.vcffile) {
        Some(r) => r,
        None => return,
    };
    let subtype_id = r.id2int("SUBTYPE");
    let sites_annotated = subtype_id >= 0 && r.info_exists(subtype_id);
    let consbp_present = r.id2int("CONSBP") >= 0;
    let mut delly_vcf = false;
    while r.read() {
        r.unpack_info();
        if !delly_vcf {
            if let Some(m) = r.get_string("SVMETHOD") {
                if m.len() >= 10 && &m[..10] == "EMBL.DELLY" && consbp_present {
                    delly_vcf = true;
                }
            }
        }
        if !delly_vcf {
            eprintln!("Error: Delly genotyping requires local SV assembly (INFO/CONSENSUS) and breakpoint (INFO/CONSBP) introduced in delly v1.1.7!");
            break;
        }
        let mut sv = SVRecord::default();
        let tid = hd.name2id(&r.rid_name());
        sv.chr = tid;
        sv.sv_start = r.pos() as i32 + 1;
        sv.id = svs.len() as i32;
        sv.mapq = r.qual() as i32;
        let ref_allele = r.allele(0);
        if !ref_allele.is_empty() && ref_allele != "." {
            let alt = r.allele(1);
            if !alt.is_empty() && alt != "." {
                sv.alleles = format!("{},{}", ref_allele, alt);
            }
        }
        match (r.get_string("SVTYPE"), r.get_string("CT")) {
            (Some(svt), Some(ct)) => sv.svt = decode_orientation(&ct, &svt),
            _ => continue,
        }
        sv.precise = r.get_flag("PRECISE");
        sv.pe_support = match r.get_i32("PE") {
            Some(v) => v[0],
            None => {
                if sv.precise {
                    0
                } else {
                    2
                }
            }
        };
        if sv.svt != 4 {
            sv.ins_len = r.get_i32("INSLEN").map(|v| v[0]).unwrap_or(0);
        } else {
            match r.get_i32("SVLEN") {
                Some(v) => sv.ins_len = v[0],
                None => continue,
            }
        }
        sv.hom_len = r.get_i32("HOMLEN").map(|v| v[0]).unwrap_or(0);
        if sites_annotated {
            sv.anno.from_sites = true;
            sv.anno.hom_len = sv.hom_len;
            if let Some(sub) = r.get_string("SUBTYPE") {
                match sub.as_str() {
                    "INS:ME:ALU" => sv.anno.seq_type = 1,
                    "INS:ME:LINE1" => sv.anno.seq_type = 2,
                    "INS:ME:SVA" => sv.anno.seq_type = 3,
                    "INS:NUMT" => sv.anno.seq_type = 4,
                    "INS:LTR" => sv.anno.seq_type = 5,
                    "INS:HERVK" => sv.anno.seq_type = 6,
                    "INS:TR" | "DEL:TR" => {
                        if let Some(p) = r.get_i32("TRPERIOD") {
                            sv.anno.seq_type = 7;
                            sv.anno.tr_period = p[0];
                            if let Some(cp) = r.get_f32("TRCOPIES") {
                                sv.anno.tr_copies = cp[0];
                            }
                        }
                    }
                    _ => {}
                }
                if sv.anno.seq_type >= 1 && sv.anno.seq_type <= 6 {
                    if let Some(s) = r.get_string("INSSTRAND") {
                        sv.anno.is_rc = s.starts_with('-');
                    }
                }
            }
        }
        sv.sr_support = r.get_i32("SR").map(|v| v[0]).unwrap_or(0);
        sv.chr2 = tid;
        sv.sv_end = r.pos() as i32 + 1;
        if sv.svt < DELLY_SVT_TRANS {
            if let Some(e) = r.get_i32("END") {
                sv.sv_end = e[0];
            }
            if !sv.alleles.is_empty() {
                let (sym, cp) = is_symbolic(&sv.alleles);
                if !sym {
                    sv.sv_end = sv.sv_start + cp.unwrap() as i32;
                }
            }
        } else {
            if let Some(c2) = r.get_string("CHR2") {
                sv.chr2 = hd.name2id(&c2);
            }
            if let Some(p2) = r.get_i32("POS2") {
                sv.sv_end = p2[0];
            }
        }
        if let Some(cons) = r.get_string("CONSENSUS") {
            sv.consensus = cons;
            if let Some(cb) = r.get_i32("CONSBP") {
                sv.cons_bp = cb[0];
            }
        } else {
            sv.precise = false;
        }
        match r.get_i32("CIPOS") {
            Some(v) => {
                sv.ciposlow = v[0];
                sv.ciposhigh = v[1];
            }
            None => {
                sv.ciposlow = -50;
                sv.ciposhigh = 50;
            }
        }
        match r.get_i32("CIEND") {
            Some(v) => {
                sv.ciendlow = v[0];
                sv.ciendhigh = v[1];
            }
            None => {
                sv.ciendlow = -50;
                sv.ciendhigh = 50;
            }
        }
        sv.pe_map_quality = r.get_i32("MAPQ").map(|v| v[0] as u8 as i32).unwrap_or(0);
        sv.sr_map_quality = r.get_i32("SRMAPQ").map(|v| v[0] as u8 as i32).unwrap_or(0);
        sv.sr_align_quality = r.get_f32("SRQ").map(|v| v[0]).unwrap_or(0.0);
        if let Some(a) = r.get_i32("ALLELEID") {
            sv.alleleid = a[0];
            if let Some(n) = r.get_i32("NALLELE") {
                sv.nallele = n[0];
            }
            if let (Some(ts), Some(te)) = (r.get_i32("TRSTART"), r.get_i32("TREND")) {
                sv.tr_start = ts[0] - 1;
                sv.tr_end = te[0];
            }
        }
        svs.push(sv);
    }
}
