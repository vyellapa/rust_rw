//! FFI wrapper around the bundled edlib library (compiled from delly's edlib.cpp).

#![allow(dead_code)]

use std::os::raw::{c_char, c_int};

pub const EDLIB_STATUS_OK: i32 = 0;

pub const EDLIB_EDOP_MATCH: u8 = 0;
pub const EDLIB_EDOP_INSERT: u8 = 1;
pub const EDLIB_EDOP_DELETE: u8 = 2;
pub const EDLIB_EDOP_MISMATCH: u8 = 3;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Mode {
    NW = 0,
    SHW = 1,
    HW = 2,
}

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Task {
    Distance = 0,
    Loc = 1,
    Path = 2,
}

#[repr(C)]
struct DlEdlibResult {
    status: c_int,
    edit_distance: c_int,
    end_location0: c_int,
    alignment_length: c_int,
    alignment: *mut u8,
}

extern "C" {
    fn dl_edlib_align(q: *const c_char, ql: c_int, t: *const c_char, tl: c_int, k: c_int, mode: c_int, task: c_int, out: *mut DlEdlibResult);
    fn dl_edlib_free(p: *mut u8);
}

/// Owned result of an edlib alignment
pub struct AlignResult {
    pub status: i32,
    pub edit_distance: i32,
    /// endLocations[0] (or -1 when no location was found)
    pub end_location: i32,
    pub alignment: Vec<u8>,
}

pub fn align(query: &[u8], target: &[u8], k: i32, mode: Mode, task: Task) -> AlignResult {
    let mut out = DlEdlibResult { status: 1, edit_distance: -1, end_location0: -1, alignment_length: 0, alignment: std::ptr::null_mut() };
    unsafe {
        dl_edlib_align(
            query.as_ptr() as *const c_char,
            query.len() as c_int,
            target.as_ptr() as *const c_char,
            target.len() as c_int,
            k,
            mode as c_int,
            task as c_int,
            &mut out,
        );
        let alignment = if out.alignment.is_null() {
            Vec::new()
        } else {
            let v = std::slice::from_raw_parts(out.alignment, out.alignment_length as usize).to_vec();
            dl_edlib_free(out.alignment);
            v
        };
        AlignResult { status: out.status, edit_distance: out.edit_distance, end_location: out.end_location0, alignment }
    }
}

/// infixStart (util.h)
pub fn infix_start(r: &AlignResult) -> u32 {
    let mut t_idx = r.end_location;
    for &op in &r.alignment {
        if op != EDLIB_EDOP_INSERT {
            t_idx -= 1;
        }
    }
    if t_idx >= 0 {
        (t_idx + 1) as u32
    } else {
        0
    }
}

/// infixEnd (util.h)
pub fn infix_end(r: &AlignResult) -> u32 {
    r.end_location as u32
}
