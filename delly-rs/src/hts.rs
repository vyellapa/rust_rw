//! Minimal FFI bindings to htslib (linked from the system / Homebrew installation)
//! plus small safe wrappers used throughout the port.

#![allow(non_camel_case_types, dead_code)]

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void};
use std::ptr;

pub type hts_pos_t = i64;

#[repr(C)]
pub struct htsFile {
    _p: [u8; 0],
}
#[repr(C)]
pub struct hts_idx_t {
    _p: [u8; 0],
}
#[repr(C)]
pub struct hts_itr_t {
    _p: [u8; 0],
}
#[repr(C)]
pub struct sam_hdr_t {
    _p: [u8; 0],
}
#[repr(C)]
pub struct bcf_hdr_t {
    _p: [u8; 0],
}
#[repr(C)]
pub struct bcf1_t {
    _p: [u8; 0],
}
#[repr(C)]
pub struct faidx_t {
    _p: [u8; 0],
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct bam1_core_t {
    pub pos: hts_pos_t,
    pub tid: i32,
    pub bin: u16,
    pub qual: u8,
    pub l_extranul: u8,
    pub flag: u16,
    pub l_qname: u16,
    pub n_cigar: u32,
    pub l_qseq: i32,
    pub mtid: i32,
    pub mpos: hts_pos_t,
    pub isize: hts_pos_t,
}

#[repr(C)]
pub struct bam1_t {
    pub core: bam1_core_t,
    pub id: u64,
    pub data: *mut u8,
    pub l_data: c_int,
    pub m_data: u32,
    pub mempolicy: u32,
}

pub const BAM_FPAIRED: u16 = 1;
pub const BAM_FPROPER_PAIR: u16 = 2;
pub const BAM_FUNMAP: u16 = 4;
pub const BAM_FMUNMAP: u16 = 8;
pub const BAM_FREVERSE: u16 = 16;
pub const BAM_FMREVERSE: u16 = 32;
pub const BAM_FREAD1: u16 = 64;
pub const BAM_FREAD2: u16 = 128;
pub const BAM_FSECONDARY: u16 = 256;
pub const BAM_FQCFAIL: u16 = 512;
pub const BAM_FDUP: u16 = 1024;
pub const BAM_FSUPPLEMENTARY: u16 = 2048;

pub const BAM_CMATCH: u32 = 0;
pub const BAM_CINS: u32 = 1;
pub const BAM_CDEL: u32 = 2;
pub const BAM_CREF_SKIP: u32 = 3;
pub const BAM_CSOFT_CLIP: u32 = 4;
pub const BAM_CHARD_CLIP: u32 = 5;
pub const BAM_CPAD: u32 = 6;
pub const BAM_CEQUAL: u32 = 7;
pub const BAM_CDIFF: u32 = 8;

// CRAM required fields (sam.h)
pub const SAM_QNAME: c_int = 0x00000001;
pub const SAM_FLAG: c_int = 0x00000002;
pub const SAM_RNAME: c_int = 0x00000004;
pub const SAM_POS: c_int = 0x00000008;
pub const SAM_MAPQ: c_int = 0x00000010;
pub const SAM_CIGAR: c_int = 0x00000020;
pub const SAM_RNEXT: c_int = 0x00000040;
pub const SAM_PNEXT: c_int = 0x00000080;
pub const SAM_TLEN: c_int = 0x00000100;
pub const SAM_SEQ: c_int = 0x00000200;

pub const BCF_DT_ID: c_int = 0;
pub const BCF_DT_CTG: c_int = 1;

pub const BCF_INT32_MISSING: i32 = i32::MIN;
pub const BCF_INT32_VECTOR_END: i32 = i32::MIN + 1;
pub const BCF_FLOAT_MISSING_BITS: u32 = 0x7F800001;

#[inline]
pub fn bcf_gt_unphased(a: i32) -> i32 {
    (a + 1) << 1
}
#[inline]
pub fn bcf_gt_phased(a: i32) -> i32 {
    ((a + 1) << 1) | 1
}
pub const BCF_GT_MISSING: i32 = 0;
#[inline]
pub fn bcf_gt_is_missing(v: i32) -> bool {
    (v >> 1) == 0
}
#[inline]
pub fn bcf_gt_allele(v: i32) -> i32 {
    (v >> 1) - 1
}

extern "C" {
    pub fn hts_open(fn_: *const c_char, mode: *const c_char) -> *mut htsFile;
    pub fn hts_close(fp: *mut htsFile) -> c_int;
    pub fn hts_version() -> *const c_char;
    pub fn hts_set_fai_filename(fp: *mut htsFile, fn_aux: *const c_char) -> c_int;
    pub fn hts_idx_destroy(idx: *mut hts_idx_t);
    pub fn hts_itr_destroy(itr: *mut hts_itr_t);
    pub fn hts_idx_get_stat(idx: *const hts_idx_t, tid: c_int, mapped: *mut u64, unmapped: *mut u64) -> c_int;

    pub fn sam_hdr_read(fp: *mut htsFile) -> *mut sam_hdr_t;
    pub fn sam_hdr_destroy(h: *mut sam_hdr_t);
    pub fn sam_index_load(fp: *mut htsFile, fn_: *const c_char) -> *mut hts_idx_t;
    pub fn bam_init1() -> *mut bam1_t;
    pub fn bam_destroy1(b: *mut bam1_t);

    pub fn fai_load(fn_: *const c_char) -> *mut faidx_t;
    pub fn fai_build(fn_: *const c_char) -> c_int;
    pub fn fai_destroy(fai: *mut faidx_t);
    pub fn faidx_fetch_seq(fai: *const faidx_t, c_name: *const c_char, p_beg_i: c_int, p_end_i: c_int, len: *mut c_int) -> *mut c_char;
    pub fn faidx_has_seq(fai: *const faidx_t, seq: *const c_char) -> c_int;

    pub fn bcf_hdr_init(mode: *const c_char) -> *mut bcf_hdr_t;
    pub fn bcf_hdr_destroy(h: *mut bcf_hdr_t);
    pub fn bcf_hdr_append(h: *mut bcf_hdr_t, line: *const c_char) -> c_int;
    pub fn bcf_hdr_add_sample(h: *mut bcf_hdr_t, sample: *const c_char) -> c_int;
    pub fn bcf_hdr_write(fp: *mut htsFile, h: *mut bcf_hdr_t) -> c_int;
    pub fn bcf_hdr_read(fp: *mut htsFile) -> *mut bcf_hdr_t;
    pub fn bcf_hdr_id2int(h: *const bcf_hdr_t, type_: c_int, id: *const c_char) -> c_int;
    pub fn bcf_init() -> *mut bcf1_t;
    pub fn bcf_destroy(v: *mut bcf1_t);
    pub fn bcf_clear(v: *mut bcf1_t);
    pub fn bcf_write(fp: *mut htsFile, h: *mut bcf_hdr_t, v: *mut bcf1_t) -> c_int;
    pub fn bcf_read(fp: *mut htsFile, h: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
    pub fn bcf_update_id(h: *const bcf_hdr_t, line: *mut bcf1_t, id: *const c_char) -> c_int;
    pub fn bcf_update_alleles_str(h: *const bcf_hdr_t, line: *mut bcf1_t, alleles: *const c_char) -> c_int;
    pub fn bcf_update_filter(h: *const bcf_hdr_t, line: *mut bcf1_t, flt_ids: *mut c_int, n: c_int) -> c_int;
    pub fn bcf_update_format_string(h: *const bcf_hdr_t, line: *mut bcf1_t, key: *const c_char, values: *mut *const c_char, n: c_int) -> c_int;
    pub fn bcf_index_build(fn_: *const c_char, min_shift: c_int) -> c_int;

    // shim.c
    pub fn dl_sam_itr_next(fp: *mut htsFile, itr: *mut hts_itr_t, r: *mut bam1_t) -> c_int;
    pub fn dl_sam_itr_queryi(idx: *const hts_idx_t, tid: c_int, beg: hts_pos_t, end: hts_pos_t) -> *mut hts_itr_t;
    pub fn dl_hdr_n_targets(h: *const sam_hdr_t) -> c_int;
    pub fn dl_hdr_target_name(h: *const sam_hdr_t, i: c_int) -> *const c_char;
    pub fn dl_hdr_target_len(h: *const sam_hdr_t, i: c_int) -> u32;
    pub fn dl_hdr_text(h: *const sam_hdr_t) -> *const c_char;
    pub fn dl_hdr_l_text(h: *const sam_hdr_t) -> usize;
    pub fn dl_bam_name2id(h: *mut sam_hdr_t, r: *const c_char) -> c_int;
    pub fn dl_hts_set_required_fields(fp: *mut htsFile, fields: c_int) -> c_int;
    pub fn dl_hts_set_threads(fp: *mut htsFile, n: c_int) -> c_int;

    pub fn dl_bcf_hdr_nsamples(h: *const bcf_hdr_t) -> c_int;
    pub fn dl_bcf_hdr_name2id(h: *const bcf_hdr_t, id: *const c_char) -> c_int;
    pub fn dl_bcf_hdr_id2name(h: *const bcf_hdr_t, rid: c_int) -> *const c_char;
    pub fn dl_bcf_hdr_idinfo_exists_info(h: *const bcf_hdr_t, id: c_int) -> c_int;
    pub fn dl_bcf_hdr_id2int_id(h: *const bcf_hdr_t, id: *const c_char) -> c_int;
    pub fn dl_bcf_set_rid(r: *mut bcf1_t, rid: c_int);
    pub fn dl_bcf_set_pos(r: *mut bcf1_t, pos: hts_pos_t);
    pub fn dl_bcf_set_qual(r: *mut bcf1_t, q: f32);
    pub fn dl_bcf_get_rid(r: *const bcf1_t) -> c_int;
    pub fn dl_bcf_get_pos(r: *const bcf1_t) -> hts_pos_t;
    pub fn dl_bcf_get_qual(r: *const bcf1_t) -> f32;
    pub fn dl_bcf_n_allele(r: *const bcf1_t) -> c_int;
    pub fn dl_bcf_allele(r: *const bcf1_t, i: c_int) -> *const c_char;
    pub fn dl_bcf_update_info_int32(h: *const bcf_hdr_t, r: *mut bcf1_t, key: *const c_char, v: *const i32, n: c_int) -> c_int;
    pub fn dl_bcf_update_info_float(h: *const bcf_hdr_t, r: *mut bcf1_t, key: *const c_char, v: *const f32, n: c_int) -> c_int;
    pub fn dl_bcf_update_info_flag(h: *const bcf_hdr_t, r: *mut bcf1_t, key: *const c_char) -> c_int;
    pub fn dl_bcf_update_info_string(h: *const bcf_hdr_t, r: *mut bcf1_t, key: *const c_char, s: *const c_char) -> c_int;
    pub fn dl_bcf_update_format_int32(h: *const bcf_hdr_t, r: *mut bcf1_t, key: *const c_char, v: *const i32, n: c_int) -> c_int;
    pub fn dl_bcf_update_genotypes(h: *const bcf_hdr_t, r: *mut bcf1_t, v: *const i32, n: c_int) -> c_int;
    pub fn dl_bcf_get_info_int32(h: *const bcf_hdr_t, r: *mut bcf1_t, tag: *const c_char, dst: *mut *mut i32, ndst: *mut c_int) -> c_int;
    pub fn dl_bcf_get_info_float(h: *const bcf_hdr_t, r: *mut bcf1_t, tag: *const c_char, dst: *mut *mut f32, ndst: *mut c_int) -> c_int;
    pub fn dl_bcf_get_info_string(h: *const bcf_hdr_t, r: *mut bcf1_t, tag: *const c_char, dst: *mut *mut c_char, ndst: *mut c_int) -> c_int;
    pub fn dl_bcf_get_info_flag(h: *const bcf_hdr_t, r: *mut bcf1_t, tag: *const c_char) -> c_int;
    pub fn dl_bcf_unpack_info(r: *mut bcf1_t) -> c_int;
}

pub fn cstr(s: &str) -> CString {
    CString::new(s).expect("string with interior NUL")
}

pub fn hts_version_string() -> String {
    unsafe { CStr::from_ptr(hts_version()).to_string_lossy().into_owned() }
}

// ---------------------------------------------------------------------------
// BAM header
// ---------------------------------------------------------------------------

/// Owned copy of the information we need from a SAM/BAM header.
#[derive(Clone, Debug)]
pub struct Header {
    pub names: Vec<String>,
    pub lens: Vec<u32>,
    pub text: String,
}

impl Header {
    pub fn n_targets(&self) -> i32 {
        self.names.len() as i32
    }
    /// bam_name2id: exact-name lookup, -1 if missing
    pub fn name2id(&self, name: &str) -> i32 {
        self.names.iter().position(|n| n == name).map(|p| p as i32).unwrap_or(-1)
    }
    unsafe fn from_raw(h: *const sam_hdr_t) -> Header {
        let n = dl_hdr_n_targets(h);
        let mut names = Vec::with_capacity(n as usize);
        let mut lens = Vec::with_capacity(n as usize);
        for i in 0..n {
            names.push(CStr::from_ptr(dl_hdr_target_name(h, i)).to_string_lossy().into_owned());
            lens.push(dl_hdr_target_len(h, i));
        }
        let tp = dl_hdr_text(h);
        let text = if tp.is_null() {
            String::new()
        } else {
            let l = dl_hdr_l_text(h);
            String::from_utf8_lossy(std::slice::from_raw_parts(tp as *const u8, l)).into_owned()
        };
        Header { names, lens, text }
    }
}

// ---------------------------------------------------------------------------
// Indexed BAM/CRAM reader
// ---------------------------------------------------------------------------

pub struct IndexedBam {
    fp: *mut htsFile,
    idx: *mut hts_idx_t,
    hdr: *mut sam_hdr_t,
    pub header: Header,
}

unsafe impl Send for IndexedBam {}

impl IndexedBam {
    /// sam_open + hts_set_fai_filename + sam_index_load + sam_hdr_read (in delly's order)
    pub fn open(path: &str, genome: Option<&str>) -> Option<IndexedBam> {
        let p = cstr(path);
        unsafe {
            let fp = hts_open(p.as_ptr(), b"r\0".as_ptr() as *const c_char);
            if fp.is_null() {
                return None;
            }
            if let Some(g) = genome {
                let g = cstr(g);
                hts_set_fai_filename(fp, g.as_ptr());
            }
            let idx = sam_index_load(fp, p.as_ptr());
            if idx.is_null() {
                hts_close(fp);
                return None;
            }
            let hdr = sam_hdr_read(fp);
            if hdr.is_null() {
                hts_idx_destroy(idx);
                hts_close(fp);
                return None;
            }
            let header = Header::from_raw(hdr);
            Some(IndexedBam { fp, idx, hdr, header })
        }
    }

    /// Open a file and only read its header (no index needed)
    pub fn read_header(path: &str) -> Option<Header> {
        let p = cstr(path);
        unsafe {
            let fp = hts_open(p.as_ptr(), b"r\0".as_ptr() as *const c_char);
            if fp.is_null() {
                return None;
            }
            let hdr = sam_hdr_read(fp);
            if hdr.is_null() {
                hts_close(fp);
                return None;
            }
            let h = Header::from_raw(hdr);
            sam_hdr_destroy(hdr);
            hts_close(fp);
            Some(h)
        }
    }

    pub fn set_threads(&mut self, n: i32) {
        if n > 1 {
            unsafe {
                dl_hts_set_threads(self.fp, n);
            }
        }
    }

    pub fn set_required_fields(&mut self, fields: c_int) {
        unsafe {
            dl_hts_set_required_fields(self.fp, fields);
        }
    }

    /// hts_idx_get_stat mapped count
    pub fn mapped(&self, tid: i32) -> u64 {
        let mut m: u64 = 0;
        let mut u: u64 = 0;
        unsafe {
            hts_idx_get_stat(self.idx, tid, &mut m, &mut u);
        }
        m
    }

    /// Iterate all records overlapping [beg, end) of tid; the callback returns false to stop early.
    pub fn for_each_in<F: FnMut(&Record) -> bool>(&mut self, tid: i32, beg: i64, end: i64, mut f: F) -> bool {
        unsafe {
            let iter = dl_sam_itr_queryi(self.idx, tid, beg, end);
            if iter.is_null() {
                return false;
            }
            let b = bam_init1();
            let rec = Record { b };
            while dl_sam_itr_next(self.fp, iter, b) >= 0 {
                if !f(&rec) {
                    break;
                }
            }
            bam_destroy1(b);
            hts_itr_destroy(iter);
        }
        true
    }
}

impl Drop for IndexedBam {
    fn drop(&mut self) {
        unsafe {
            sam_hdr_destroy(self.hdr);
            hts_idx_destroy(self.idx);
            hts_close(self.fp);
        }
    }
}

/// Borrowed view over a bam1_t
pub struct Record {
    b: *mut bam1_t,
}

impl Record {
    #[inline]
    fn core(&self) -> &bam1_core_t {
        unsafe { &(*self.b).core }
    }
    #[inline]
    pub fn pos(&self) -> i64 {
        self.core().pos
    }
    #[inline]
    pub fn tid(&self) -> i32 {
        self.core().tid
    }
    #[inline]
    pub fn mtid(&self) -> i32 {
        self.core().mtid
    }
    #[inline]
    pub fn mpos(&self) -> i64 {
        self.core().mpos
    }
    #[inline]
    pub fn isize(&self) -> i64 {
        self.core().isize
    }
    #[inline]
    pub fn qual(&self) -> u8 {
        self.core().qual
    }
    #[inline]
    pub fn flag(&self) -> u16 {
        self.core().flag
    }
    #[inline]
    pub fn l_qseq(&self) -> i32 {
        self.core().l_qseq
    }
    #[inline]
    pub fn n_cigar(&self) -> u32 {
        self.core().n_cigar
    }
    /// Read name without the trailing NUL
    #[inline]
    pub fn qname(&self) -> &[u8] {
        unsafe {
            let d = (*self.b).data;
            let l = self.core().l_qname as usize;
            let s = std::slice::from_raw_parts(d, l);
            let n = s.iter().position(|&c| c == 0).unwrap_or(l);
            &s[..n]
        }
    }
    /// Raw CIGAR words
    #[inline]
    pub fn cigar(&self) -> &[u32] {
        unsafe {
            let d = (*self.b).data.add(self.core().l_qname as usize) as *const u32;
            std::slice::from_raw_parts(d, self.core().n_cigar as usize)
        }
    }
    /// Packed 4-bit sequence
    #[inline]
    fn seq_ptr(&self) -> *const u8 {
        unsafe { (*self.b).data.add(self.core().l_qname as usize + (self.core().n_cigar as usize) * 4) }
    }
    /// Decoded sequence using "=ACMGRSVTWYHKDBN"
    pub fn seq_string(&self) -> Vec<u8> {
        const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let n = self.l_qseq().max(0) as usize;
        let mut out = Vec::with_capacity(n);
        let p = self.seq_ptr();
        unsafe {
            for i in 0..n {
                let byte = *p.add(i >> 1);
                let code = (byte >> ((!i & 1) << 2)) & 0xf;
                out.push(NT16[code as usize]);
            }
        }
        out
    }
}

#[inline]
pub fn cigar_op(c: u32) -> u32 {
    c & 0xf
}
#[inline]
pub fn cigar_oplen(c: u32) -> u32 {
    c >> 4
}

// ---------------------------------------------------------------------------
// FASTA index
// ---------------------------------------------------------------------------

pub struct Faidx {
    fai: *mut faidx_t,
}
unsafe impl Send for Faidx {}
unsafe impl Sync for Faidx {}

impl Faidx {
    pub fn load(path: &str) -> Option<Faidx> {
        let p = cstr(path);
        let fai = unsafe { fai_load(p.as_ptr()) };
        if fai.is_null() {
            None
        } else {
            Some(Faidx { fai })
        }
    }
    pub fn build(path: &str) -> bool {
        let p = cstr(path);
        unsafe { fai_build(p.as_ptr()) != -1 }
    }
    pub fn has_seq(&self, name: &str) -> bool {
        let n = cstr(name);
        unsafe { faidx_has_seq(self.fai, n.as_ptr()) != 0 }
    }
    /// faidx_fetch_seq(name, 0, len) → whole sequence.  The returned buffer keeps a
    /// trailing NUL byte (as in C) so code that reads one past the end behaves alike.
    pub fn fetch(&self, name: &str, beg: i32, end: i32) -> Option<RefSeq> {
        let n = cstr(name);
        let mut len: c_int = -1;
        unsafe {
            let s = faidx_fetch_seq(self.fai, n.as_ptr(), beg, end, &mut len);
            if s.is_null() {
                return None;
            }
            let l = if len > 0 { len as usize } else { 0 };
            let mut v = Vec::with_capacity(l + 1);
            v.extend_from_slice(std::slice::from_raw_parts(s as *const u8, l));
            v.push(0);
            libc::free(s as *mut c_void);
            Some(RefSeq { data: v })
        }
    }
}

impl Drop for Faidx {
    fn drop(&mut self) {
        unsafe { fai_destroy(self.fai) }
    }
}

/// Reference sequence with a trailing NUL sentinel
pub struct RefSeq {
    data: Vec<u8>,
}

impl RefSeq {
    /// Sequence without the sentinel
    #[inline]
    pub fn seq(&self) -> &[u8] {
        &self.data[..self.data.len() - 1]
    }
    /// Byte access including the NUL sentinel at index len
    #[inline]
    pub fn at(&self, i: usize) -> u8 {
        self.data[i]
    }
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len() - 1
    }
}

// ---------------------------------------------------------------------------
// BCF writer
// ---------------------------------------------------------------------------

pub struct BcfWriter {
    pub fp: *mut htsFile,
    pub hdr: *mut bcf_hdr_t,
    pub rec: *mut bcf1_t,
}

impl BcfWriter {
    pub fn open(path: &str, mode: &str) -> Option<BcfWriter> {
        let p = cstr(path);
        let m = cstr(mode);
        unsafe {
            let fp = hts_open(p.as_ptr(), m.as_ptr());
            if fp.is_null() {
                return None;
            }
            let hdr = bcf_hdr_init(b"w\0".as_ptr() as *const c_char);
            Some(BcfWriter { fp, hdr, rec: ptr::null_mut() })
        }
    }
    pub fn hdr_append(&mut self, line: &str) {
        let l = cstr(line);
        unsafe {
            bcf_hdr_append(self.hdr, l.as_ptr());
        }
    }
    pub fn add_sample(&mut self, s: Option<&str>) {
        unsafe {
            match s {
                Some(s) => {
                    let c = cstr(s);
                    bcf_hdr_add_sample(self.hdr, c.as_ptr());
                }
                None => {
                    bcf_hdr_add_sample(self.hdr, ptr::null());
                }
            }
        }
    }
    pub fn write_header(&mut self) -> bool {
        unsafe { bcf_hdr_write(self.fp, self.hdr) == 0 }
    }
    pub fn nsamples(&self) -> i32 {
        unsafe { dl_bcf_hdr_nsamples(self.hdr) }
    }
    pub fn id2int(&self, id: &str) -> i32 {
        let c = cstr(id);
        unsafe { bcf_hdr_id2int(self.hdr, BCF_DT_ID, c.as_ptr()) }
    }
    pub fn name2id(&self, id: &str) -> i32 {
        let c = cstr(id);
        unsafe { dl_bcf_hdr_name2id(self.hdr, c.as_ptr()) }
    }
    pub fn init_rec(&mut self) {
        self.rec = unsafe { bcf_init() };
    }
    pub fn set_rid(&mut self, rid: i32) {
        unsafe { dl_bcf_set_rid(self.rec, rid) }
    }
    pub fn set_pos(&mut self, pos: i64) {
        unsafe { dl_bcf_set_pos(self.rec, pos) }
    }
    pub fn set_qual(&mut self, q: f32) {
        unsafe { dl_bcf_set_qual(self.rec, q) }
    }
    pub fn update_id(&mut self, id: &str) {
        let c = cstr(id);
        unsafe {
            bcf_update_id(self.hdr, self.rec, c.as_ptr());
        }
    }
    pub fn update_alleles_str(&mut self, a: &str) {
        let c = cstr(a);
        unsafe {
            bcf_update_alleles_str(self.hdr, self.rec, c.as_ptr());
        }
    }
    pub fn update_filter(&mut self, id: i32) {
        let mut v = id;
        unsafe {
            bcf_update_filter(self.hdr, self.rec, &mut v, 1);
        }
    }
    pub fn info_flag(&mut self, key: &str) {
        let k = cstr(key);
        unsafe {
            dl_bcf_update_info_flag(self.hdr, self.rec, k.as_ptr());
        }
    }
    pub fn info_string(&mut self, key: &str, val: &str) {
        let k = cstr(key);
        let v = cstr(val);
        unsafe {
            dl_bcf_update_info_string(self.hdr, self.rec, k.as_ptr(), v.as_ptr());
        }
    }
    pub fn info_i32(&mut self, key: &str, vals: &[i32]) {
        let k = cstr(key);
        unsafe {
            dl_bcf_update_info_int32(self.hdr, self.rec, k.as_ptr(), vals.as_ptr(), vals.len() as c_int);
        }
    }
    pub fn info_f32(&mut self, key: &str, vals: &[f32]) {
        let k = cstr(key);
        unsafe {
            dl_bcf_update_info_float(self.hdr, self.rec, k.as_ptr(), vals.as_ptr(), vals.len() as c_int);
        }
    }
    pub fn format_i32(&mut self, key: &str, vals: &[i32]) {
        let k = cstr(key);
        unsafe {
            dl_bcf_update_format_int32(self.hdr, self.rec, k.as_ptr(), vals.as_ptr(), vals.len() as c_int);
        }
    }
    pub fn genotypes(&mut self, vals: &[i32]) {
        unsafe {
            dl_bcf_update_genotypes(self.hdr, self.rec, vals.as_ptr(), vals.len() as c_int);
        }
    }
    pub fn format_string(&mut self, key: &str, vals: &[String]) {
        let k = cstr(key);
        let cs: Vec<CString> = vals.iter().map(|s| cstr(s)).collect();
        let mut ps: Vec<*const c_char> = cs.iter().map(|c| c.as_ptr()).collect();
        unsafe {
            bcf_update_format_string(self.hdr, self.rec, k.as_ptr(), ps.as_mut_ptr(), ps.len() as c_int);
        }
    }
    pub fn write_rec(&mut self) {
        unsafe {
            bcf_write(self.fp, self.hdr, self.rec);
            bcf_clear(self.rec);
        }
    }
    pub fn close(mut self) {
        unsafe {
            if !self.rec.is_null() {
                bcf_destroy(self.rec);
                self.rec = ptr::null_mut();
            }
            bcf_hdr_destroy(self.hdr);
            hts_close(self.fp);
        }
        std::mem::forget(self);
    }
}

pub fn bcf_index_build_csi(path: &str, min_shift: i32) {
    let p = cstr(path);
    unsafe {
        bcf_index_build(p.as_ptr(), min_shift);
    }
}

// ---------------------------------------------------------------------------
// BCF reader (for genotyping mode)
// ---------------------------------------------------------------------------

pub struct BcfReader {
    fp: *mut htsFile,
    pub hdr: *mut bcf_hdr_t,
    pub rec: *mut bcf1_t,
    ibuf: *mut i32,
    ni: c_int,
    fbuf: *mut f32,
    nf: c_int,
    sbuf: *mut c_char,
    ns: c_int,
}

impl BcfReader {
    pub fn open(path: &str) -> Option<BcfReader> {
        let p = cstr(path);
        unsafe {
            let fp = hts_open(p.as_ptr(), b"r\0".as_ptr() as *const c_char);
            if fp.is_null() {
                return None;
            }
            let hdr = bcf_hdr_read(fp);
            if hdr.is_null() {
                hts_close(fp);
                return None;
            }
            Some(BcfReader {
                fp,
                hdr,
                rec: bcf_init(),
                ibuf: ptr::null_mut(),
                ni: 0,
                fbuf: ptr::null_mut(),
                nf: 0,
                sbuf: ptr::null_mut(),
                ns: 0,
            })
        }
    }
    pub fn read(&mut self) -> bool {
        unsafe { bcf_read(self.fp, self.hdr, self.rec) == 0 }
    }
    pub fn unpack_info(&mut self) {
        unsafe {
            dl_bcf_unpack_info(self.rec);
        }
    }
    pub fn id2int(&self, id: &str) -> i32 {
        let c = cstr(id);
        unsafe { dl_bcf_hdr_id2int_id(self.hdr, c.as_ptr()) }
    }
    pub fn info_exists(&self, id: i32) -> bool {
        unsafe { dl_bcf_hdr_idinfo_exists_info(self.hdr, id) != 0 }
    }
    pub fn rid_name(&self) -> String {
        unsafe { CStr::from_ptr(dl_bcf_hdr_id2name(self.hdr, dl_bcf_get_rid(self.rec))).to_string_lossy().into_owned() }
    }
    pub fn pos(&self) -> i64 {
        unsafe { dl_bcf_get_pos(self.rec) }
    }
    pub fn qual(&self) -> f32 {
        unsafe { dl_bcf_get_qual(self.rec) }
    }
    pub fn allele(&self, i: i32) -> String {
        unsafe { CStr::from_ptr(dl_bcf_allele(self.rec, i)).to_string_lossy().into_owned() }
    }
    pub fn n_allele(&self) -> i32 {
        unsafe { dl_bcf_n_allele(self.rec) }
    }
    pub fn get_i32(&mut self, tag: &str) -> Option<Vec<i32>> {
        let t = cstr(tag);
        unsafe {
            let r = dl_bcf_get_info_int32(self.hdr, self.rec, t.as_ptr(), &mut self.ibuf, &mut self.ni);
            if r > 0 {
                Some(std::slice::from_raw_parts(self.ibuf, r as usize).to_vec())
            } else {
                None
            }
        }
    }
    pub fn get_f32(&mut self, tag: &str) -> Option<Vec<f32>> {
        let t = cstr(tag);
        unsafe {
            let r = dl_bcf_get_info_float(self.hdr, self.rec, t.as_ptr(), &mut self.fbuf, &mut self.nf);
            if r > 0 {
                Some(std::slice::from_raw_parts(self.fbuf, r as usize).to_vec())
            } else {
                None
            }
        }
    }
    pub fn get_string(&mut self, tag: &str) -> Option<String> {
        let t = cstr(tag);
        unsafe {
            let r = dl_bcf_get_info_string(self.hdr, self.rec, t.as_ptr(), &mut self.sbuf, &mut self.ns);
            if r > 0 {
                Some(CStr::from_ptr(self.sbuf).to_string_lossy().into_owned())
            } else {
                None
            }
        }
    }
    pub fn get_flag(&mut self, tag: &str) -> bool {
        let t = cstr(tag);
        unsafe { dl_bcf_get_info_flag(self.hdr, self.rec, t.as_ptr()) > 0 }
    }
}

impl Drop for BcfReader {
    fn drop(&mut self) {
        unsafe {
            libc::free(self.ibuf as *mut c_void);
            libc::free(self.fbuf as *mut c_void);
            libc::free(self.sbuf as *mut c_void);
            bcf_hdr_destroy(self.hdr);
            hts_close(self.fp);
            bcf_destroy(self.rec);
        }
    }
}
