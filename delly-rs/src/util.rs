//! Port of util.h helpers used by `delly sr`.

use crate::config::Config;
use crate::hts::*;
use crate::tags::*;
use std::collections::BTreeSet;
use std::io::{BufRead, Read};

/// Right-open genomic interval [lower, upper)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Interval {
    pub lower: u32,
    pub upper: u32,
}

/// boost::icl::interval_set<uint32_t> (joining) - kept sorted, disjoint, non-touching
#[derive(Clone, Debug, Default)]
pub struct IntervalSet {
    pub ivs: Vec<Interval>,
}

impl IntervalSet {
    pub fn is_empty(&self) -> bool {
        self.ivs.is_empty()
    }
    /// Insert a right-open interval, joining overlapping or touching intervals
    pub fn insert(&mut self, lo: u32, hi: u32) {
        if lo >= hi {
            return;
        }
        let mut nlo = lo;
        let mut nhi = hi;
        let mut out = Vec::with_capacity(self.ivs.len() + 1);
        let mut placed = false;
        for iv in &self.ivs {
            if iv.upper < nlo {
                out.push(*iv);
            } else if nhi < iv.lower {
                if !placed {
                    out.push(Interval { lower: nlo, upper: nhi });
                    placed = true;
                }
                out.push(*iv);
            } else {
                nlo = nlo.min(iv.lower);
                nhi = nhi.max(iv.upper);
            }
        }
        if !placed {
            out.push(Interval { lower: nlo, upper: nhi });
        }
        self.ivs = out;
    }
}

pub type ValidRegions = Vec<IntervalSet>;

/// readLength
pub fn read_length(rec: &Record) -> i32 {
    let mut hc: i32 = 0;
    for &c in rec.cigar() {
        if cigar_op(c) == BAM_CHARD_CLIP {
            hc += cigar_oplen(c) as i32;
        }
    }
    rec.l_qseq() + hc
}

/// alignmentLength
pub fn alignment_length(rec: &Record) -> u32 {
    let mut alen: u32 = 0;
    for &c in rec.cigar() {
        let op = cigar_op(c);
        if op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP {
            alen = alen.wrapping_add(cigar_oplen(c));
        }
    }
    alen
}

#[inline]
pub fn half_alignment_length(rec: &Record) -> u32 {
    alignment_length(rec) / 2
}

/// boost::to_upper_copy for ASCII
#[inline]
pub fn to_upper(s: &[u8]) -> Vec<u8> {
    s.iter().map(|c| c.to_ascii_uppercase()).collect()
}

/// reverseComplement (in place); uppercases, leaves unknown characters as the original
/// character at that (reversed) position, exactly as util.h does.
pub fn reverse_complement(seq: &mut [u8]) {
    let rev: Vec<u8> = seq.iter().rev().map(|c| c.to_ascii_uppercase()).collect();
    for (i, &c) in rev.iter().enumerate() {
        match c {
            b'A' => seq[i] = b'T',
            b'C' => seq[i] = b'G',
            b'G' => seq[i] = b'C',
            b'T' => seq[i] = b'A',
            b'N' => seq[i] = b'N',
            _ => {}
        }
    }
}

/// entropy (clang fuses `ent += freq * log2(freq)` into an FMA)
pub fn entropy(st: &[u8]) -> f64 {
    let alphabet: BTreeSet<i8> = st.iter().map(|&c| c as i8).collect();
    let mut ent: f64 = 0.0;
    for &c in &alphabet {
        let ctr = st.iter().filter(|&&s| s as i8 == c).count() as i32;
        let freq = ctr as f64 / st.len() as f64;
        ent = freq.mul_add(freq.log2(), ent);
    }
    -ent
}

/// getSMTag
pub fn get_sm_tag(header: &str, file_name: &str) -> String {
    let mut sm: BTreeSet<String> = BTreeSet::new();
    let mut rg_present = false;
    for line in header.split('\n') {
        if line.starts_with("@RG") {
            for kv in line.split('\t') {
                if let Some(sp) = kv.find(':') {
                    if &kv[..sp] == "SM" {
                        rg_present = true;
                        sm.insert(kv[sp + 1..].to_string());
                    }
                }
            }
        }
    }
    if !rg_present {
        file_name.to_string()
    } else if sm.len() == 1 {
        sm.iter().next().unwrap().clone()
    } else if sm.len() > 1 {
        eprintln!("Warning: Multiple sample names (@RG:SM) present in the BAM file!");
        sm.iter().next().unwrap().clone()
    } else {
        "unknown".to_string()
    }
}

/// checkSampleNames
pub fn check_sample_names(c: &mut Config) {
    let mut ucount: u32 = 0;
    let mut snames: BTreeSet<String> = BTreeSet::new();
    for i in 0..c.files.len() {
        while snames.contains(&c.sample_name[i]) {
            eprintln!("Warning: Duplicate sample names: {}", c.sample_name[i]);
            c.sample_name[i] = format!("{}_{}", c.sample_name[i], ucount);
            ucount += 1;
            eprintln!("Warning: Changing sample name to {}", c.sample_name[i]);
        }
        snames.insert(c.sample_name[i].clone());
    }
}

/// Parse an integer like boost::lexical_cast<int32_t> (whole token must be a valid int)
fn lexical_cast_i32(s: &str) -> Option<i32> {
    if s.is_empty() {
        return None;
    }
    let b = s.as_bytes();
    let (neg, digits) = match b[0] {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let v: i64 = digits.parse::<i64>().ok()?;
    let v = if neg { -v } else { v };
    if v < i32::MIN as i64 || v > i32::MAX as i64 {
        None
    } else {
        Some(v as i32)
    }
}

fn is_gz(path: &str) -> bool {
    let mut b = [0u8; 2];
    match std::fs::File::open(path) {
        Ok(mut f) => f.read_exact(&mut b).is_ok() && b[0] == 0x1f && b[1] == 0x8b,
        Err(_) => false,
    }
}

/// _parseExcludeIntervals
pub fn parse_exclude_intervals(c: &Config, hdr: &Header, valid: &mut ValidRegions) -> bool {
    let n = hdr.n_targets() as usize;
    valid.clear();
    valid.resize(n, IntervalSet::default());
    let mut exclg: Vec<IntervalSet> = vec![IntervalSet::default(); n];
    let mut valid_chr = vec![true; n];
    if c.has_exclude_file {
        let lines: Vec<String> = if is_gz(&c.exclude) {
            match gunzip_lines(&c.exclude) {
                Some(l) => l,
                None => Vec::new(),
            }
        } else {
            match std::fs::File::open(&c.exclude) {
                Ok(f) => std::io::BufReader::new(f).split(b'\n').map(|l| String::from_utf8_lossy(&l.unwrap_or_default()).into_owned()).collect(),
                Err(_) => Vec::new(),
            }
        };
        for line in lines {
            let toks: Vec<&str> = line.split(|ch| ch == ' ' || ch == '\t' || ch == ',' || ch == ';').filter(|t| !t.is_empty()).collect();
            if toks.is_empty() {
                continue;
            }
            let tid = hdr.name2id(toks[0]);
            if tid >= 0 {
                if toks.len() > 1 {
                    let start = match lexical_cast_i32(toks[1]) {
                        Some(v) => v,
                        None => {
                            eprintln!("Exclude file needs to be in tab-delimited format: chr, start, end");
                            eprintln!("Offending line: {}", line);
                            return false;
                        }
                    };
                    if toks.len() > 2 {
                        let end = match lexical_cast_i32(toks[2]) {
                            Some(v) => v,
                            None => {
                                eprintln!("Exclude file needs to be in tab-delimited format: chr, start, end");
                                eprintln!("Offending line: {}", line);
                                return false;
                            }
                        };
                        if start < end {
                            exclg[tid as usize].insert(start as u32, end as u32);
                        } else {
                            eprintln!("Exclude file needs to be in tab-delimited format (chr, start, end) and start < end.");
                            eprintln!("Offending line: {}", line);
                            return false;
                        }
                    } else {
                        eprintln!("Exclude file needs to be in tab-delimited format: chr, start, end");
                        eprintln!("Offending line: {}", line);
                        return false;
                    }
                } else {
                    valid_chr[tid as usize] = false;
                }
            }
        }
    }
    for i in 0..n {
        if !valid_chr[i] {
            continue;
        }
        let mut istart: u32 = 0;
        for iv in &exclg[i].ivs {
            if istart.wrapping_add(1) < iv.lower {
                valid[i].insert(istart, iv.lower - 1);
            }
            istart = iv.upper;
        }
        if istart.wrapping_add(1) < hdr.lens[i] {
            valid[i].insert(istart, hdr.lens[i]);
        }
    }
    true
}

/// Minimal gzip line reader for exclude files (single or multi-member gzip via flate in zlib).
fn gunzip_lines(path: &str) -> Option<Vec<String>> {
    // Use htslib's bgzf/hts reader which transparently handles gzip
    let p = cstr(path);
    unsafe {
        extern "C" {
            fn hts_open(fn_: *const std::os::raw::c_char, mode: *const std::os::raw::c_char) -> *mut htsFile;
            fn hts_getline(fp: *mut htsFile, delimiter: std::os::raw::c_int, str_: *mut KString) -> std::os::raw::c_int;
            fn hts_close(fp: *mut htsFile) -> std::os::raw::c_int;
        }
        #[repr(C)]
        struct KString {
            l: usize,
            m: usize,
            s: *mut std::os::raw::c_char,
        }
        let fp = hts_open(p.as_ptr(), b"r\0".as_ptr() as *const std::os::raw::c_char);
        if fp.is_null() {
            return None;
        }
        let mut ks = KString { l: 0, m: 0, s: std::ptr::null_mut() };
        let mut out = Vec::new();
        while hts_getline(fp, b'\n' as i32, &mut ks) >= 0 {
            let sl = std::slice::from_raw_parts(ks.s as *const u8, ks.l);
            out.push(String::from_utf8_lossy(sl).into_owned());
        }
        if !ks.s.is_null() {
            libc::free(ks.s as *mut std::os::raw::c_void);
        }
        hts_close(fp);
        Some(out)
    }
}

/// getVariability
pub fn get_variability(lib: &[LibraryInfo]) -> i32 {
    let mut v = 0;
    for l in lib {
        if l.max_normal_isize > v {
            v = l.max_normal_isize;
        }
        if l.rs > v {
            v = l.rs;
        }
    }
    v
}

/// medianVector: nth_element then v[n]; the value is uniquely defined
pub fn median_vector(v: &mut [u8]) -> i32 {
    let n = v.len() / 2;
    v.select_nth_unstable(n);
    v[n] as i32
}

/// getLibraryParams for one file
pub fn library_params(c: &Config, valid: &ValidRegions, file_c: usize) -> LibraryInfo {
    let mut lib = LibraryInfo::default();
    let mut bam = match IndexedBam::open(&c.files[file_c], Some(&c.genome)) {
        Some(b) => b,
        None => return lib,
    };
    let max_alignments_screened: u32 = 10_000_000;
    let max_num_alignments: u32 = 1_000_000;
    let min_num_alignments: u32 = 1000;
    let mut alignment_count: u32 = 0;
    let mut processed_num_pairs: u32 = 0;
    let mut processed_num_reads: u32 = 0;
    let mut rplus: u32 = 0;
    let mut nonrplus: u32 = 0;
    let mut vec_isize: Vec<u32> = Vec::with_capacity(min_num_alignments as usize);
    let mut read_size: Vec<u32> = Vec::with_capacity(min_num_alignments as usize);
    let mut lib_characterized = false;
    let n = bam.header.n_targets();
    'outer: for ref_index in 0..n {
        if valid[ref_index as usize].is_empty() {
            continue;
        }
        let ivs = valid[ref_index as usize].ivs.clone();
        for iv in ivs {
            if lib_characterized {
                break;
            }
            bam.for_each_in(ref_index, iv.lower as i64, iv.upper as i64, |rec| {
                if (rec.flag() & BAM_FREAD2 == 0) && (rec.l_qseq() < 65000) {
                    if rec.flag() & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP) != 0 {
                        return true;
                    }
                    if (alignment_count > max_alignments_screened)
                        || ((processed_num_reads >= max_num_alignments) && (processed_num_pairs == 0))
                        || (processed_num_pairs >= max_num_alignments)
                    {
                        lib_characterized = true;
                        return false;
                    }
                    alignment_count += 1;
                    if processed_num_reads < max_num_alignments {
                        read_size.push(rec.l_qseq() as u32);
                        processed_num_reads += 1;
                    }
                    if (rec.flag() & BAM_FPAIRED != 0) && (rec.flag() & BAM_FMUNMAP == 0) && (rec.tid() == rec.mtid()) && processed_num_pairs < max_num_alignments {
                        // abs(hts_pos_t) then narrowing into vector<uint32_t>
                        vec_isize.push(rec.isize().wrapping_abs() as u32);
                        if get_sv_type(rec) == 2 {
                            rplus += 1;
                        } else {
                            nonrplus += 1;
                        }
                        processed_num_pairs += 1;
                    }
                }
                true
            });
            if lib_characterized {
                break 'outer;
            }
        }
        if lib_characterized {
            break;
        }
    }
    if processed_num_reads >= min_num_alignments {
        read_size.sort_unstable();
        lib.rs = read_size[read_size.len() / 2] as i32;
    }
    if processed_num_pairs >= min_num_alignments {
        vec_isize.sort_unstable();
        let median = vec_isize[vec_isize.len() / 2] as i32;
        let mut abs_dev: Vec<u32> = vec_isize.iter().map(|&v| ((v as i32) - median).unsigned_abs()).collect();
        abs_dev.sort_unstable();
        let mad = abs_dev[abs_dev.len() / 2] as i32;
        if (50..=100000).contains(&median) {
            if rplus < nonrplus {
                eprintln!("Warning: Sample has a non-default paired-end layout! File: {}", c.files[file_c]);
                eprintln!("The expected paired-end orientation is   ---Read1--->      <---Read2---  which is the default illumina paired-end layout.");
            } else {
                lib.median = median;
                lib.mad = mad;
                lib.max_normal_isize = median + (c.mad_normal_cutoff as i32 * mad);
                lib.min_normal_isize = median - (c.mad_normal_cutoff as i32 * mad);
                if lib.min_normal_isize < 0 {
                    lib.min_normal_isize = 0;
                }
                lib.max_isize_cutoff = median + (c.mad_cutoff as i32 * mad);
                lib.min_isize_cutoff = median - (c.mad_cutoff as i32 * mad);
                lib.max_isize_cutoff = lib.max_isize_cutoff.max(2 * lib.rs);
                lib.max_isize_cutoff = lib.max_isize_cutoff.max(500);
                if lib.min_isize_cutoff < 0 {
                    lib.min_isize_cutoff = 0;
                }
            }
        }
    }
    lib
}

/// Timestamp prefix like boost::posix_time::to_simple_string(second_clock::local_time())
pub fn now_string() -> String {
    unsafe {
        let t = libc::time(std::ptr::null_mut());
        let mut tm: libc::tm = std::mem::zeroed();
        libc::localtime_r(&t, &mut tm);
        const MON: [&str; 12] = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
        format!("{:04}-{}-{:02} {:02}:{:02}:{:02}", tm.tm_year + 1900, MON[tm.tm_mon as usize], tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec)
    }
}

/// boost::gregorian::to_iso_string(today) → YYYYMMDD
pub fn today_iso() -> String {
    unsafe {
        let t = libc::time(std::ptr::null_mut());
        let mut tm: libc::tm = std::mem::zeroed();
        libc::localtime_r(&t, &mut tm);
        format!("{:04}{:02}{:02}", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday)
    }
}

static START: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();

pub fn log(msg: &str) {
    eprintln!("[{}] {}", now_string(), msg);
    if std::env::var_os("DELLY_TIMING").is_some() {
        let t = START.get_or_init(std::time::Instant::now);
        eprintln!("  [timing] {:.3}s", t.elapsed().as_secs_f64());
    }
}
