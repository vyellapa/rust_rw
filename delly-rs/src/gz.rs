//! gzip writer matching boost::iostreams::gzip_compressor (default params) using the
//! system zlib, as linked by the reference binary.

use std::os::raw::{c_char, c_int, c_uint, c_ulong, c_void};
use std::io::{self, Write};

#[repr(C)]
struct ZStream {
    next_in: *const u8,
    avail_in: c_uint,
    total_in: c_ulong,
    next_out: *mut u8,
    avail_out: c_uint,
    total_out: c_ulong,
    msg: *const c_char,
    state: *mut c_void,
    zalloc: *const c_void,
    zfree: *const c_void,
    opaque: *mut c_void,
    data_type: c_int,
    adler: c_ulong,
    reserved: c_ulong,
}

#[link(name = "z")]
extern "C" {
    fn deflateInit2_(strm: *mut ZStream, level: c_int, method: c_int, window_bits: c_int, mem_level: c_int, strategy: c_int, version: *const c_char, stream_size: c_int) -> c_int;
    fn deflate(strm: *mut ZStream, flush: c_int) -> c_int;
    fn deflateEnd(strm: *mut ZStream) -> c_int;
    fn crc32(crc: c_ulong, buf: *const u8, len: c_uint) -> c_ulong;
    fn zlibVersion() -> *const c_char;
}

const Z_NO_FLUSH: c_int = 0;
const Z_OK: c_int = 0;
const Z_FINISH: c_int = 4;
const Z_STREAM_END: c_int = 1;

/// Compress `data` into gzip format like boost's gzip_compressor and write it to `path`.
pub fn write_gzip_file(path: &str, data: &[u8]) -> io::Result<()> {
    let mut out = io::BufWriter::new(std::fs::File::create(path)?);
    write_gzip(&mut out, data, 1 << 20)?;
    out.flush()
}

fn write_gzip(out: &mut impl Write, data: &[u8], chunk_size: usize) -> io::Result<()> {
    assert!(chunk_size > 0);
    unsafe {
        let mut s: ZStream = std::mem::zeroed();
        // level -1 (Z_DEFAULT_COMPRESSION), Z_DEFLATED, raw deflate (-15), memLevel 8, default strategy
        let init = deflateInit2_(&mut s, -1, 8, -15, 8, 0, zlibVersion(), std::mem::size_of::<ZStream>() as c_int);
        if init != Z_OK {
            return Err(io::Error::other(format!("zlib initialization failed: {}", init)));
        }
        // Keep the initialized stream at this address until deflateEnd, including on I/O errors.
        let result = (|| {
            out.write_all(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff])?;
            let mut buf = vec![0u8; 1 << 16];
            let mut crc = 0;
            for chunk in data.chunks(chunk_size.min(c_uint::MAX as usize)) {
                crc = crc32(crc, chunk.as_ptr(), chunk.len() as c_uint);
                s.next_in = chunk.as_ptr();
                s.avail_in = chunk.len() as c_uint;
                while s.avail_in > 0 {
                    s.next_out = buf.as_mut_ptr();
                    s.avail_out = buf.len() as c_uint;
                    let r = deflate(&mut s, Z_NO_FLUSH);
                    if r != Z_OK {
                        return Err(io::Error::other(format!("zlib compression failed: {}", r)));
                    }
                    out.write_all(&buf[..buf.len() - s.avail_out as usize])?;
                }
            }
            loop {
                s.next_out = buf.as_mut_ptr();
                s.avail_out = buf.len() as c_uint;
                let r = deflate(&mut s, Z_FINISH);
                if r != Z_OK && r != Z_STREAM_END {
                    return Err(io::Error::other(format!("zlib finish failed: {}", r)));
                }
                out.write_all(&buf[..buf.len() - s.avail_out as usize])?;
                if r == Z_STREAM_END {
                    break;
                }
            }
            out.write_all(&(crc as u32).to_le_bytes())?;
            // gzip ISIZE deliberately stores the total length modulo 2^32.
            out.write_all(&(data.len() as u32).to_le_bytes())
        })();
        deflateEnd(&mut s);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunking_does_not_change_gzip_bytes() {
        let mut seed = 1u64;
        let data: Vec<u8> = (0..200_000).map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 32) as u8
        }).collect();
        for input in [&[][..], &b"hello hello hello"[..], &data[..]] {
            let mut expected = Vec::new();
            write_gzip(&mut expected, input, usize::MAX).unwrap();
            for size in [1, 127, 65536] {
                let mut got = Vec::new();
                write_gzip(&mut got, input, size).unwrap();
                assert_eq!(got, expected, "input chunk size {}", size);
            }
        }
    }

    #[test]
    fn write_failures_are_returned() {
        struct FailedWriter;
        impl Write for FailedWriter {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> { Err(io::Error::other("full")) }
            fn flush(&mut self) -> io::Result<()> { Ok(()) }
        }
        assert!(write_gzip(&mut FailedWriter, b"test", 2).is_err());
    }

    #[cfg(all(unix, target_pointer_width = "64"))]
    #[test]
    #[ignore = "compresses and decompresses over 4 GiB; requires gzip"]
    fn input_over_four_gib_is_not_truncated() {
        use std::process::{Command, Stdio};
        let len = (1usize << 32) + 3;
        unsafe {
            let p = libc::mmap(std::ptr::null_mut(), len, libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_PRIVATE | libc::MAP_ANON, -1, 0);
            assert_ne!(p, libc::MAP_FAILED);
            let data = std::slice::from_raw_parts_mut(p as *mut u8, len);
            data[..3].copy_from_slice(b"ABC");
            data[len - 3..].copy_from_slice(b"XYZ");
            let path = std::env::temp_dir().join(format!("delly-large-gzip-{}.gz", std::process::id()));
            let result = write_gzip_file(path.to_str().unwrap(), data);
            libc::munmap(p, len);
            result.unwrap();
            let mut child = Command::new("gzip").arg("-dc").arg(&path).stdout(Stdio::piped()).spawn().unwrap();
            let count = io::copy(&mut child.stdout.take().unwrap(), &mut io::sink()).unwrap();
            let status = child.wait().unwrap();
            std::fs::remove_file(path).unwrap();
            assert!(status.success(), "gzip CRC/length validation failed");
            assert_eq!(count, len as u64);
        }
    }
}
