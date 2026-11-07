//! Configuration and command-line parsing for `delly sr` (boost::program_options look-alike).

use crate::ploidy::SexModel;
use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug)]
pub struct DnaScore {
    pub match_: i32,
    pub mismatch: i32,
    pub go: i32,
    pub ge: i32,
    pub inf: i32,
}

impl DnaScore {
    pub fn new(m: i32, mm: i32, go: i32, ge: i32) -> Self {
        DnaScore { match_: m, mismatch: mm, go, ge, inf: 1_000_000 }
    }
}

impl Default for DnaScore {
    fn default() -> Self {
        DnaScore::new(5, -4, -10, -1)
    }
}

#[derive(Clone, Debug)]
pub struct Config {
    pub min_map_qual: u16,
    pub min_tra_qual: u16,
    pub min_geno_qual: u16,
    pub mad_cutoff: u16,
    pub mad_normal_cutoff: u16,
    pub nchr: i32,
    pub minimum_flank_size: i32,
    pub indelsize: i32,
    pub geno_cap: i32,
    pub min_cons_window: i32,
    pub graph_pruning: u32,
    pub min_ref_sep: u32,
    pub max_read_sep: u32,
    pub min_clip: u32,
    pub max_geno_read_count: u32,
    pub min_clique_size: u32,
    pub max_threads: u32,
    pub max_read_per_sv: u32,
    pub flank_quality: f32,
    pub mei_min_frac: f32,
    pub tr_min_frac: f32,
    pub has_exclude_file: bool,
    pub has_vcf_file: bool,
    pub has_dump_file: bool,
    pub svtset: BTreeSet<i32>,
    pub aliscore: DnaScore,
    pub outfile: String,
    pub vcffile: String,
    pub genome: String,
    pub exclude: String,
    pub dumpfile: String,
    pub files: Vec<String>,
    pub sample_name: Vec<String>,
    pub sex_arg: String,
    pub sex_model: SexModel,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            min_map_qual: 1,
            min_tra_qual: 20,
            min_geno_qual: 5,
            mad_cutoff: 9,
            mad_normal_cutoff: 5,
            nchr: 0,
            minimum_flank_size: 13,
            indelsize: 1000,
            geno_cap: 25,
            min_cons_window: 100,
            graph_pruning: 1000,
            min_ref_sep: 25,
            max_read_sep: 40,
            min_clip: 25,
            max_geno_read_count: 250,
            min_clique_size: 2,
            max_threads: 4,
            max_read_per_sv: 20,
            flank_quality: 0.95,
            mei_min_frac: 0.8,
            tr_min_frac: 0.85,
            has_exclude_file: false,
            has_vcf_file: false,
            has_dump_file: false,
            svtset: BTreeSet::new(),
            aliscore: DnaScore::default(),
            outfile: String::new(),
            vcffile: String::new(),
            genome: String::new(),
            exclude: String::new(),
            dumpfile: String::new(),
            files: Vec::new(),
            sample_name: Vec::new(),
            sex_arg: "auto".to_string(),
            sex_model: SexModel::default(),
        }
    }
}

/// Option specification
struct Opt {
    long: &'static str,
    short: Option<char>,
    takes_value: bool,
    visible: bool,
    desc: &'static str,
    default: &'static str,
}

const OPTS: &[Opt] = &[
    Opt { long: "help", short: Some('?'), takes_value: false, visible: true, desc: "show help message", default: "" },
    Opt { long: "svtype", short: Some('t'), takes_value: true, visible: true, desc: "SV type to compute [DEL, INS, DUP, INV, BND, ALL]", default: "ALL" },
    Opt { long: "genome", short: Some('g'), takes_value: true, visible: true, desc: "genome fasta file", default: "" },
    Opt { long: "exclude", short: Some('x'), takes_value: true, visible: true, desc: "file with regions to exclude", default: "" },
    Opt { long: "outfile", short: Some('o'), takes_value: true, visible: true, desc: "BCF output file", default: "" },
    Opt { long: "threads", short: Some('h'), takes_value: true, visible: true, desc: "number of threads", default: "4" },
    Opt { long: "map-qual", short: Some('q'), takes_value: true, visible: true, desc: "min. paired-end (PE) mapping quality", default: "1" },
    Opt { long: "qual-tra", short: Some('r'), takes_value: true, visible: true, desc: "min. PE quality for translocation", default: "20" },
    Opt { long: "mad-cutoff", short: Some('s'), takes_value: true, visible: true, desc: "insert size cutoff, median+s*MAD (deletions only)", default: "9" },
    Opt { long: "minclip", short: Some('c'), takes_value: true, visible: true, desc: "min. clipping length", default: "25" },
    Opt { long: "min-clique-size", short: Some('z'), takes_value: true, visible: true, desc: "min. PE/SR clique size", default: "2" },
    Opt { long: "minrefsep", short: Some('m'), takes_value: true, visible: true, desc: "min. reference separation", default: "25" },
    Opt { long: "maxreadsep", short: Some('n'), takes_value: true, visible: true, desc: "max. read separation", default: "40" },
    Opt { long: "max-reads", short: Some('p'), takes_value: true, visible: true, desc: "max. reads for consensus computation", default: "20" },
    Opt { long: "vcffile", short: Some('v'), takes_value: true, visible: true, desc: "input VCF/BCF file for genotyping", default: "" },
    Opt { long: "geno-qual", short: Some('u'), takes_value: true, visible: true, desc: "min. mapping quality for genotyping", default: "5" },
    Opt { long: "dump", short: Some('d'), takes_value: true, visible: true, desc: "gzipped output file for SV-reads (optional)", default: "" },
    Opt { long: "max-geno-count", short: Some('a'), takes_value: true, visible: true, desc: "max. number of reads aligned for SR genotyping", default: "250" },
    Opt { long: "sex", short: None, takes_value: true, visible: true, desc: "sample sex [auto, male, female, none, file]", default: "auto" },
    Opt { long: "input-file", short: None, takes_value: true, visible: false, desc: "input file", default: "" },
    Opt { long: "pruning", short: Some('j'), takes_value: true, visible: false, desc: "PE graph pruning cutoff", default: "1000" },
    Opt { long: "cons-window", short: Some('w'), takes_value: true, visible: false, desc: "consensus window", default: "100" },
];

pub struct Parsed {
    pub values: std::collections::BTreeMap<&'static str, Vec<String>>,
}

impl Parsed {
    pub fn count(&self, k: &str) -> bool {
        self.values.contains_key(k)
    }
    pub fn get(&self, k: &str) -> Option<&str> {
        self.values.get(k).and_then(|v| v.last()).map(|s| s.as_str())
    }
}

fn find_long(name: &str) -> Result<&'static Opt, String> {
    if let Some(o) = OPTS.iter().find(|o| o.long == name) {
        return Ok(o);
    }
    // boost allow_guessing: unique prefix
    let cands: Vec<&Opt> = OPTS.iter().filter(|o| o.long.starts_with(name)).collect();
    if cands.len() == 1 {
        Ok(cands[0])
    } else if cands.is_empty() {
        Err(format!("unrecognised option '--{}'", name))
    } else {
        Err(format!("option '--{}' is ambiguous", name))
    }
}

pub fn parse_args(args: &[String]) -> Result<Parsed, String> {
    let mut values: std::collections::BTreeMap<&'static str, Vec<String>> = std::collections::BTreeMap::new();
    let mut i = 1;
    while i < args.len() {
        let a = &args[i];
        if a == "--" {
            for p in &args[i + 1..] {
                values.entry("input-file").or_default().push(p.clone());
            }
            break;
        } else if let Some(rest) = a.strip_prefix("--") {
            let (name, inline) = match rest.find('=') {
                Some(p) => (&rest[..p], Some(rest[p + 1..].to_string())),
                None => (rest, None),
            };
            let o = find_long(name)?;
            if o.takes_value {
                let v = match inline {
                    Some(v) => v,
                    None => {
                        i += 1;
                        if i >= args.len() {
                            return Err(format!("the required argument for option '--{}' is missing", o.long));
                        }
                        args[i].clone()
                    }
                };
                values.entry(o.long).or_default().push(v);
            } else {
                values.entry(o.long).or_default().push(String::new());
            }
        } else if a.len() > 1 && a.starts_with('-') {
            let chars: Vec<char> = a[1..].chars().collect();
            let mut k = 0;
            while k < chars.len() {
                let ch = chars[k];
                let o = OPTS.iter().find(|o| o.short == Some(ch)).ok_or_else(|| format!("unrecognised option '-{}'", ch))?;
                if o.takes_value {
                    let rest: String = chars[k + 1..].iter().collect();
                    let v = if !rest.is_empty() {
                        rest
                    } else {
                        i += 1;
                        if i >= args.len() {
                            return Err(format!("the required argument for option '--{}' is missing", o.long));
                        }
                        args[i].clone()
                    };
                    values.entry(o.long).or_default().push(v);
                    break;
                } else {
                    values.entry(o.long).or_default().push(String::new());
                    k += 1;
                }
            }
        } else {
            values.entry("input-file").or_default().push(a.clone());
        }
        i += 1;
    }
    Ok(Parsed { values })
}

pub fn print_usage(cmd: &str) {
    eprintln!();
    eprintln!("Usage: delly {} [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ...", cmd);
    let groups: [(&str, &[&str]); 3] = [
        ("Generic options", &["help", "svtype", "genome", "exclude", "outfile", "threads"]),
        ("Discovery options", &["map-qual", "qual-tra", "mad-cutoff", "minclip", "min-clique-size", "minrefsep", "maxreadsep", "max-reads"]),
        ("Genotyping options", &["vcffile", "geno-qual", "dump", "max-geno-count", "sex"]),
    ];
    for (title, names) in groups.iter() {
        eprintln!("{}:", title);
        for n in names.iter() {
            let o = OPTS.iter().find(|o| o.long == *n).unwrap();
            if !o.visible {
                continue;
            }
            let mut left = match o.short {
                Some(s) => format!("  -{} [ --{} ]", s, o.long),
                None => format!("  --{}", o.long),
            };
            if o.takes_value {
                if o.default.is_empty() {
                    left.push_str(" arg");
                } else {
                    left.push_str(&format!(" arg (={})", o.default));
                }
            }
            eprintln!("{:<39} {}", left, o.desc);
        }
        eprintln!();
    }
}

fn num<T: std::str::FromStr>(p: &Parsed, k: &str, default: T) -> Result<T, String> {
    match p.get(k) {
        Some(v) => v.trim().parse::<T>().map_err(|_| format!("the argument ('{}') for option '--{}' is invalid", v, k)),
        None => Ok(default),
    }
}

/// Fill the config from parsed options; returns (config, svtype string)
pub fn config_from(p: &Parsed) -> Result<(Config, String), String> {
    let mut c = Config::default();
    let svtype = p.get("svtype").unwrap_or("ALL").to_string();
    c.genome = p.get("genome").unwrap_or("").to_string();
    c.exclude = p.get("exclude").unwrap_or("").to_string();
    c.outfile = p.get("outfile").unwrap_or("").to_string();
    c.max_threads = num(p, "threads", 4u32)?;
    c.min_map_qual = num(p, "map-qual", 1u16)?;
    c.min_tra_qual = num(p, "qual-tra", 20u16)?;
    c.mad_cutoff = num(p, "mad-cutoff", 9u16)?;
    c.min_clip = num(p, "minclip", 25u32)?;
    c.min_clique_size = num(p, "min-clique-size", 2u32)?;
    c.min_ref_sep = num(p, "minrefsep", 25u32)?;
    c.max_read_sep = num(p, "maxreadsep", 40u32)?;
    c.max_read_per_sv = num(p, "max-reads", 20u32)?;
    c.vcffile = p.get("vcffile").unwrap_or("").to_string();
    c.min_geno_qual = num(p, "geno-qual", 5u16)?;
    c.dumpfile = p.get("dump").unwrap_or("").to_string();
    c.max_geno_read_count = num(p, "max-geno-count", 250u32)?;
    c.sex_arg = p.get("sex").unwrap_or("auto").to_string();
    c.files = p.values.get("input-file").cloned().unwrap_or_default();
    c.graph_pruning = num(p, "pruning", 1000u32)?;
    c.min_cons_window = num(p, "cons-window", 100i32)?;
    Ok((c, svtype))
}

/// _svTypesToCompute
pub fn sv_types_to_compute(c: &mut Config, svtype: &str) -> bool {
    use crate::tags::DELLY_SVT_TRANS as T;
    if svtype == "ALL" {
        return true;
    }
    for tok in svtype.split(',').filter(|t| !t.is_empty()) {
        match tok {
            "DEL" => {
                c.svtset.insert(2);
            }
            "INS" => {
                c.svtset.insert(4);
            }
            "DUP" => {
                c.svtset.insert(3);
            }
            "INV" => {
                c.svtset.insert(0);
                c.svtset.insert(1);
            }
            "INV_3to3" => {
                c.svtset.insert(0);
            }
            "INV_5to5" => {
                c.svtset.insert(1);
            }
            "BND" => {
                for k in 0..4 {
                    c.svtset.insert(T + k);
                }
            }
            "BND_3to3" => {
                c.svtset.insert(T);
            }
            "BND_5to5" => {
                c.svtset.insert(T + 1);
            }
            "BND_3to5" => {
                c.svtset.insert(T + 2);
            }
            "BND_5to3" => {
                c.svtset.insert(T + 3);
            }
            _ => return false,
        }
    }
    true
}

impl Config {
    /// true if svt should be computed (empty set = all)
    #[inline]
    pub fn want(&self, svt: i32) -> bool {
        self.svtset.is_empty() || self.svtset.contains(&svt)
    }
}
