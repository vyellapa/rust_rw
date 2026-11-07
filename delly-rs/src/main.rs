//! delly-rs: a parallel Rust port of Delly's short-read SV caller (`delly sr`).
//!
//! The BCF output is byte-identical to Delly v2.8.1 built with Apple clang/libc++,
//! Boost 1.92 and the same htslib.

mod align;
mod assemble;
mod cluster;
mod config;
mod coverage;
mod cxx_hash;
mod cxx_sort;
mod edlib;
mod gz;
mod hts;
mod junction;
mod mei;
mod msa;
mod ploidy;
mod scan;
mod split;
mod svanno;
mod tags;
mod util;
mod vcf;

use config::*;
use cxx_sort::cpp_sort;
use hts::*;
use rayon::prelude::*;
use std::path::Path;
use tags::*;
use util::*;

const VERSION: &str = "2.8.1";

fn print_title(title: &str) {
    eprintln!("**********************************************************************");
    eprintln!("Program: Delly");
    eprintln!("This is free software, and you are welcome to redistribute it under");
    eprintln!("certain conditions (BSD License); for license details use '-l'.");
    eprintln!("This program comes with ABSOLUTELY NO WARRANTY; for details use '-w'.");
    eprintln!();
    eprintln!("{} (Version: {})", title, VERSION);
    eprintln!("Contact: Tobias Rausch (rausch@embl.de)");
    eprintln!("**********************************************************************");
    eprintln!();
}

fn display_usage() {
    eprintln!("Usage: delly <command> <arguments>");
    eprintln!();
    eprintln!("SV calling:");
    eprintln!("    sr           short-read SV discovery and genotyping");
    eprintln!();
}

fn file_ok(p: &str) -> bool {
    match std::fs::metadata(p) {
        Ok(m) => m.is_file() && m.len() > 0,
        Err(_) => false,
    }
}

/// boost::filesystem::path::stem()
fn path_stem(p: &str) -> String {
    let name = Path::new(p).file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    if name == "." || name == ".." {
        return name;
    }
    match name.rfind('.') {
        Some(0) | None => name,
        Some(i) => name[..i].to_string(),
    }
}

/// _outfileValid
fn outfile_valid(outfile: &str) -> bool {
    let p = Path::new(outfile);
    let outdir = match p.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
        _ => std::env::current_dir().unwrap_or_default(),
    };
    if !outdir.exists() {
        eprintln!("Output directory does not exist: \"{}\"", outdir.display());
        return false;
    }
    if std::fs::File::create(outfile).is_err() || !Path::new(outfile).is_file() {
        eprintln!("Fail to open output file {}", outfile);
        return false;
    }
    let _ = std::fs::remove_file(outfile);
    true
}

fn delly_run(c: &mut Config) -> i32 {
    let hdr = match IndexedBam::read_header(&c.files[0]) {
        Some(h) => h,
        None => {
            eprintln!("Fail to open header for {}", c.files[0]);
            return 1;
        }
    };
    let mut valid: ValidRegions = Vec::new();
    if !parse_exclude_intervals(c, &hdr, &mut valid) {
        eprintln!("Delly couldn't parse exclude intervals!");
        return 1;
    }

    // Library parameters (independent per sample)
    let cr: &Config = c;
    let mut sample_lib: Vec<LibraryInfo> = (0..cr.files.len()).into_par_iter().map(|f| library_params(cr, &valid, f)).collect();
    for (f, l) in sample_lib.iter().enumerate() {
        if l.rs == 0 {
            eprintln!("Sample has not enough data to estimate library parameters! File: {}", c.files[f]);
            return 1;
        }
    }

    let mut svs: Vec<SVRecord>;
    if !c.has_vcf_file {
        let (mut pe, mut sr, sr_store) = scan::scan_pe_and_sr(c, &valid, &hdr, &mut sample_lib);
        assemble::assemble_split_reads(c, &valid, &hdr, &sr_store, &mut sr);
        drop(sr_store);
        scan::merge_sort(&mut pe, &mut sr);
        svs = pe;
    } else {
        svs = Vec::new();
        vcf::vcf_parse(c, &hdr, &mut svs);
    }

    // Re-number SVs
    cpp_sort(&mut svs);
    for (i, sv) in svs.iter_mut().enumerate() {
        sv.id = i as i32;
    }

    // Genotyping
    let ann = if !svs.is_empty() {
        match coverage::annotate_coverage(c, &sample_lib, &mut svs, &hdr) {
            Ok(ann) => ann,
            Err(e) => {
                eprintln!("Fail to write dump file {}: {}", c.dumpfile, e);
                return 1;
            }
        }
    } else { coverage::Annotation { rc: Vec::new(), jct: Vec::new(), span: Vec::new() } };
    vcf::vcf_output(c, &hdr, &svs, &ann);

    log("Sample statistics");
    for f in 0..c.files.len() {
        let l = &sample_lib[f];
        eprintln!("Sample:{},ReadSize={},Median={},MAD={},UniqueDiscordantPairs={}", c.sample_name[f], l.rs, l.median, l.mad, l.abnormal_pairs);
    }
    log("Done.");
    0
}

fn delly(args: &[String]) -> i32 {
    let parsed = match parse_args(args) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };
    let (mut c, svtype) = match config_from(&parsed) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };
    if parsed.count("help") || !parsed.count("input-file") || !parsed.count("genome") {
        print_usage(&args[0]);
        return 0;
    }

    let mut cmd = String::from("delly ");
    for a in args {
        cmd.push_str(a);
        cmd.push(' ');
    }
    eprintln!("[{}] {}", now_string(), cmd);

    if !sv_types_to_compute(&mut c, &svtype) {
        eprintln!("Please specify a valid SV type, i.e., -t INV or -t DEL,INV without spaces.");
        return 1;
    }
    let hw = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) as u32;
    c.max_threads = c.max_threads.max(1).min(hw);
    c.has_dump_file = parsed.count("dump");
    if c.min_clique_size < 2 {
        c.min_clique_size = 2;
    }
    if c.min_map_qual > c.min_tra_qual {
        c.min_tra_qual = c.min_map_qual;
    }

    if !file_ok(&c.genome) {
        eprintln!("Reference file is missing: {}", c.genome);
        return 1;
    }
    if Faidx::load(&c.genome).is_none() {
        if !Faidx::build(&c.genome) {
            eprintln!("Fail to open genome fai index for {}", c.genome);
            return 1;
        }
    }

    // Check input files
    let nf = c.files.len();
    c.sample_name = vec![String::new(); nf];
    let mut inferred_sex = vec![0u8; nf];
    let mut x_ratio = vec![-1.0f64; nf];
    let mut y_ratio = vec![-1.0f64; nf];
    c.nchr = 0;
    let fai = Faidx::load(&c.genome).expect("genome index");
    for f in 0..nf {
        if !file_ok(&c.files[f]) {
            eprintln!("Alignment file is missing: {}", c.files[f]);
            return 1;
        }
        let bam = match IndexedBam::open(&c.files[f], None) {
            Some(b) => b,
            None => {
                eprintln!("Fail to open index for {}", c.files[f]);
                return 1;
            }
        };
        let h = bam.header.clone();
        drop(bam);
        if c.nchr == 0 {
            c.nchr = h.n_targets();
        } else if c.nchr != h.n_targets() {
            eprintln!("BAM files have different number of chromosomes!");
            return 1;
        }
        for name in &h.names {
            if !fai.has_seq(name) {
                eprintln!("Warning: BAM file chromosome {} is NOT present in your reference file {} and will be skipped.", name, c.genome);
            }
        }
        c.sample_name[f] = get_sm_tag(&h.text, &path_stem(&c.files[f]));
        if f == 0 {
            ploidy::sex_chromosomes(&h, &mut c.sex_model);
        }
        if c.sex_arg == "auto" {
            let (s, xr, yr) = ploidy::infer_sex_from_reads(&c.files[f], &c.genome, &c.sex_model);
            inferred_sex[f] = s;
            x_ratio[f] = xr;
            y_ratio[f] = yr;
        }
    }
    drop(fai);
    check_sample_names(&mut c);

    if parsed.count("exclude") {
        if !file_ok(&c.exclude) {
            eprintln!("Exclude file is missing: {}", c.exclude);
            return 1;
        }
        c.has_exclude_file = true;
    } else {
        c.has_exclude_file = false;
    }

    if parsed.count("vcffile") {
        if !file_ok(&c.vcffile) {
            eprintln!("Input VCF/BCF file is missing: {}", c.vcffile);
            return 1;
        }
        if BcfReader::open(&c.vcffile).is_none() {
            eprintln!("Fail to open file {}", c.vcffile);
            return 1;
        }
        c.has_vcf_file = true;
    } else {
        c.has_vcf_file = false;
    }

    if !parsed.count("outfile") {
        c.outfile = "-".to_string();
    } else if c.outfile != "-" && !outfile_valid(&c.outfile) {
        return 1;
    }

    let names = c.sample_name.clone();
    let sex_arg = c.sex_arg.clone();
    let auto_sex = match ploidy::parse_sex(&names, &sex_arg, &mut c.sex_model) {
        Some(a) => a,
        None => return 1,
    };
    if auto_sex {
        c.sex_model.sex = inferred_sex;
    }
    if c.sex_arg != "none" && c.sex_model.x_tid != -1 {
        for f in 0..nf {
            let mut line = format!("Sample {}: sex={}", c.sample_name[f], ploidy::sex_name(c.sex_model.sex[f]));
            if x_ratio[f] >= 0.0 {
                line.push_str(&format!(" (chrX/autosome read depth {:.2}, chrY/autosome {:.2})", x_ratio[f], if y_ratio[f] < 0.0 { 0.0 } else { y_ratio[f] }));
            }
            eprintln!("{}", line);
        }
    }

    if c.min_geno_qual < 5 {
        c.min_geno_qual = 5;
    }
    c.aliscore = DnaScore::new(5, -4, -10, -1);
    c.flank_quality = 0.95;
    c.mei_min_frac = 0.8;
    c.tr_min_frac = 0.85;
    c.minimum_flank_size = 13;
    c.indelsize = 1000;
    c.geno_cap = 25;

    // Worker pool sized by --threads
    let pool = rayon::ThreadPoolBuilder::new().num_threads(c.max_threads as usize).build().expect("thread pool");
    pool.install(|| delly_run(&mut c))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        print_title("Delly");
        display_usage();
        return;
    }
    let code = match args[1].as_str() {
        "version" | "--version" | "--version-only" | "-v" => {
            eprintln!("Delly version: v{}", VERSION);
            eprintln!(" using HTSlib: v{}", hts_version_string());
            0
        }
        "help" | "--help" | "-h" | "-?" => {
            print_title("Delly");
            display_usage();
            0
        }
        "sr" => delly(&args[1..]),
        other => {
            eprintln!("Unrecognized command {}", other);
            1
        }
    };
    std::process::exit(code);
}
