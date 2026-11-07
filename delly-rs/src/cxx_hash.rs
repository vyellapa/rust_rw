//! Bit-exact Rust emulations of C++ hashing behaviour that delly's output
//! depends on.
//!
//! Reference toolchain: Apple clang 21 / libc++ (arm64 macOS, `_LIBCPP_ABI_VERSION 1`),
//! `-std=c++17 -O3 -DNDEBUG`, Boost 1.92.
//!
//! * [`libcxx_string_hash`] == `std::hash<std::string>{}(s)`.
//!   On Apple platforms `_LIBCPP_AVAILABILITY_HAS_HASH_MEMORY` is 0, so the inline
//!   header implementation `__murmur2_or_cityhash<size_t, 64>` is used (CityHash64
//!   variant), and ABI v1 does NOT define `_LIBCPP_ABI_FIX_CITYHASH_IMPLEMENTATION`
//!   (so the 4..=8 byte case truncates `a << 3` to 32 bits).
//! * [`LibcxxStringSet`] == `std::unordered_set<std::string>` built only by
//!   `insert()` on a default-constructed set: same node list, same bucket
//!   "previous node" array, same growth / rehash (`__next_prime`) and therefore the
//!   same iteration order and `bucket_count()`.
//! * Boost 1.92 `hash_range` (char, random-access), `hash_combine`, `hash_value` for
//!   integral types, and delly's `hash_string` / `hash_sr` / `hash_pair` /
//!   `hash_pair_mate`.
//!
//! Verified bit-for-bit against C++ compiled with the reference toolchain.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// libc++ std::hash<std::string>  (__functional/hash.h, __murmur2_or_cityhash<_, 64>)
// ---------------------------------------------------------------------------

const K0: u64 = 0xc3a5c85c97cb3127;
const K1: u64 = 0xb492b66fbe98f273;
const K2: u64 = 0x9ae16a3b2f90404f;
const K3: u64 = 0xc949d7c7509e6557;

#[inline]
fn load64(s: &[u8], i: usize) -> u64 {
    u64::from_le_bytes(s[i..i + 8].try_into().unwrap())
}

#[inline]
fn load32(s: &[u8], i: usize) -> u32 {
    u32::from_le_bytes(s[i..i + 4].try_into().unwrap())
}

#[inline]
fn rotate(val: u64, shift: u32) -> u64 {
    if shift == 0 { val } else { val.rotate_right(shift) }
}

#[inline]
fn shift_mix(val: u64) -> u64 {
    val ^ (val >> 47)
}

#[inline]
fn hash_len_16(u: u64, v: u64) -> u64 {
    const MUL: u64 = 0x9ddfea08eb382d69;
    let mut a = (u ^ v).wrapping_mul(MUL);
    a ^= a >> 47;
    let mut b = (v ^ a).wrapping_mul(MUL);
    b ^= b >> 47;
    b.wrapping_mul(MUL)
}

fn hash_len_0_to_16(s: &[u8]) -> u64 {
    let len = s.len() as u64;
    let n = s.len();
    if n > 8 {
        let a = load64(s, 0);
        let b = load64(s, n - 8);
        // __rotate_by_at_least_1(b + len, len)
        return hash_len_16(a, b.wrapping_add(len).rotate_right(n as u32)) ^ b;
    }
    if n >= 4 {
        let a = load32(s, 0);
        let b = load32(s, n - 4);
        // ABI v1 (no _LIBCPP_ABI_FIX_CITYHASH_IMPLEMENTATION): `__a << 3` is a
        // uint32_t shift (truncated), then added to the size_t length.
        return hash_len_16(len.wrapping_add((a << 3) as u64), b as u64);
    }
    if n > 0 {
        let a = s[0] as u32;
        let b = s[n >> 1] as u32;
        let c = s[n - 1] as u32;
        let y: u32 = a.wrapping_add(b << 8);
        let z: u32 = (len as u32).wrapping_add(c << 2);
        return shift_mix((y as u64).wrapping_mul(K2) ^ (z as u64).wrapping_mul(K3)).wrapping_mul(K2);
    }
    K2
}

fn hash_len_17_to_32(s: &[u8]) -> u64 {
    let n = s.len();
    let len = n as u64;
    let a = load64(s, 0).wrapping_mul(K1);
    let b = load64(s, 8);
    let c = load64(s, n - 8).wrapping_mul(K2);
    let d = load64(s, n - 16).wrapping_mul(K0);
    hash_len_16(
        rotate(a.wrapping_sub(b), 43)
            .wrapping_add(rotate(c, 30))
            .wrapping_add(d),
        a.wrapping_add(rotate(b ^ K3, 20))
            .wrapping_sub(c)
            .wrapping_add(len),
    )
}

#[inline]
fn weak_hash_len_32_with_seeds6(w: u64, x: u64, y: u64, z: u64, mut a: u64, mut b: u64) -> (u64, u64) {
    a = a.wrapping_add(w);
    b = rotate(b.wrapping_add(a).wrapping_add(z), 21);
    let c = a;
    a = a.wrapping_add(x);
    a = a.wrapping_add(y);
    b = b.wrapping_add(rotate(a, 44));
    (a.wrapping_add(z), b.wrapping_add(c))
}

#[inline]
fn weak_hash_len_32_with_seeds(s: &[u8], off: usize, a: u64, b: u64) -> (u64, u64) {
    weak_hash_len_32_with_seeds6(
        load64(s, off),
        load64(s, off + 8),
        load64(s, off + 16),
        load64(s, off + 24),
        a,
        b,
    )
}

fn hash_len_33_to_64(s: &[u8]) -> u64 {
    let n = s.len();
    let len = n as u64;
    let mut z = load64(s, 24);
    let mut a = load64(s, 0).wrapping_add(len.wrapping_add(load64(s, n - 16)).wrapping_mul(K0));
    let mut b = rotate(a.wrapping_add(z), 52);
    let mut c = rotate(a, 37);
    a = a.wrapping_add(load64(s, 8));
    c = c.wrapping_add(rotate(a, 7));
    a = a.wrapping_add(load64(s, 16));
    let vf = a.wrapping_add(z);
    let vs = b.wrapping_add(rotate(a, 31)).wrapping_add(c);
    a = load64(s, 16).wrapping_add(load64(s, n - 32));
    z = z.wrapping_add(load64(s, n - 8));
    b = rotate(a.wrapping_add(z), 52);
    c = rotate(a, 37);
    a = a.wrapping_add(load64(s, n - 24));
    c = c.wrapping_add(rotate(a, 7));
    a = a.wrapping_add(load64(s, n - 16));
    let wf = a.wrapping_add(z);
    let ws = b.wrapping_add(rotate(a, 31)).wrapping_add(c);
    let r = shift_mix(
        vf.wrapping_add(ws)
            .wrapping_mul(K2)
            .wrapping_add(wf.wrapping_add(vs).wrapping_mul(K0)),
    );
    shift_mix(r.wrapping_mul(K0).wrapping_add(vs)).wrapping_mul(K2)
}

/// Exactly `std::hash<std::string>{}(s)` (also `std::hash<std::string_view>`) on
/// arm64 macOS libc++ (Apple clang 21, ABI v1).
pub fn libcxx_string_hash(s: &[u8]) -> u64 {
    let n = s.len();
    if n <= 32 {
        if n <= 16 {
            return hash_len_0_to_16(s);
        } else {
            return hash_len_17_to_32(s);
        }
    } else if n <= 64 {
        return hash_len_33_to_64(s);
    }

    let mut len = n as u64;
    let mut x = load64(s, n - 40);
    let mut y = load64(s, n - 16).wrapping_add(load64(s, n - 56));
    let mut z = hash_len_16(load64(s, n - 48).wrapping_add(len), load64(s, n - 24));
    let mut v = weak_hash_len_32_with_seeds(s, n - 64, len, z);
    let mut w = weak_hash_len_32_with_seeds(s, n - 32, y.wrapping_add(K1), x);
    x = x.wrapping_mul(K1).wrapping_add(load64(s, 0));

    len = (len - 1) & !63u64;
    let mut p = 0usize;
    loop {
        x = rotate(
            x.wrapping_add(y).wrapping_add(v.0).wrapping_add(load64(s, p + 8)),
            37,
        )
        .wrapping_mul(K1);
        y = rotate(y.wrapping_add(v.1).wrapping_add(load64(s, p + 48)), 42).wrapping_mul(K1);
        x ^= w.1;
        y = y.wrapping_add(v.0).wrapping_add(load64(s, p + 40));
        z = rotate(z.wrapping_add(w.0), 33).wrapping_mul(K1);
        v = weak_hash_len_32_with_seeds(s, p, v.1.wrapping_mul(K1), x.wrapping_add(w.0));
        w = weak_hash_len_32_with_seeds(
            s,
            p + 32,
            z.wrapping_add(w.1),
            y.wrapping_add(load64(s, p + 16)),
        );
        std::mem::swap(&mut z, &mut x);
        p += 64;
        len -= 64;
        if len == 0 {
            break;
        }
    }
    hash_len_16(
        hash_len_16(v.0, w.0)
            .wrapping_add(shift_mix(y).wrapping_mul(K1))
            .wrapping_add(z),
        hash_len_16(v.1, w.1).wrapping_add(x),
    )
}

// ---------------------------------------------------------------------------
// libc++ __next_prime (src/hash.cpp) and __hash_table helpers
// ---------------------------------------------------------------------------

const SMALL_PRIMES: [u64; 48] = [
    0, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
    193, 197, 199, 211,
];

// Residues mod 210 that are coprime to 2*3*5*7.
const INDICES: [u64; 48] = [
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181,
    187, 191, 193, 197, 199, 209,
];

/// `std::__next_prime(n)` from libc++ `hash.cpp`: the smallest prime `>= n`
/// (with `__next_prime(0) == 0`, `__next_prime(1) == 2`).
pub fn libcxx_next_prime(n: u64) -> u64 {
    const L: u64 = 210;
    if n <= SMALL_PRIMES[SMALL_PRIMES.len() - 1] {
        // *lower_bound(small_primes, n)
        let i = SMALL_PRIMES.partition_point(|&p| p < n);
        return SMALL_PRIMES[i];
    }
    assert!(n <= 0xFFFF_FFFF_FFFF_FFC5, "__next_prime overflow");
    let mut k0 = n / L;
    let mut idx = INDICES.partition_point(|&r| r < n - k0 * L);
    let mut cand = L * k0 + INDICES[idx];
    loop {
        if is_prime_wheel(cand) {
            return cand;
        }
        idx += 1;
        if idx == INDICES.len() {
            k0 += 1;
            idx = 0;
        }
        cand = L * k0 + INDICES[idx];
    }
}

/// Primality for a candidate known to be coprime to 210 and > 211
/// (same answer as libc++'s trial division loop).
fn is_prime_wheel(n: u64) -> bool {
    // Divide by 11..=199 (small primes, j = 5 .. N-2), then by numbers coprime to 210.
    for &p in &SMALL_PRIMES[5..SMALL_PRIMES.len() - 1] {
        let q = n / p;
        if q < p {
            return true;
        }
        if n == q * p {
            return false;
        }
    }
    let mut base = 210u64;
    loop {
        for &r in &INDICES {
            let i = base + r;
            let q = n / i;
            if q < i {
                return true;
            }
            if n == q * i {
                return false;
            }
        }
        base += 210;
    }
}

#[inline]
fn is_hash_power2(bc: usize) -> bool {
    bc > 2 && (bc & (bc - 1)) == 0
}

#[inline]
fn constrain_hash(h: u64, bc: usize) -> usize {
    let bc64 = bc as u64;
    if bc & bc.wrapping_sub(1) == 0 {
        (h & bc64.wrapping_sub(1)) as usize
    } else if h < bc64 {
        h as usize
    } else {
        (h % bc64) as usize
    }
}

#[inline]
fn next_hash_pow2(n: usize) -> usize {
    if n < 2 { n } else { 1usize << (usize::BITS - (n - 1).leading_zeros()) }
}

// ---------------------------------------------------------------------------
// std::unordered_set emulation (unique keys, insert-only)
// ---------------------------------------------------------------------------

const NIL: usize = usize::MAX;
const BEFORE_BEGIN: usize = 0; // index of libc++'s `__first_node_` (p1)

struct Node<K> {
    next: usize,
    hash: u64,
    value: Option<K>, // None only for the before-begin sentinel
}

/// Generic emulation of libc++ `std::unordered_set<K>` (unique keys) whose hasher is
/// `std::hash<std::string>` over the key bytes. Built only via `insert`.
pub struct LibcxxHashSet<K: AsRef<[u8]>> {
    nodes: Vec<Node<K>>, // nodes[0] = __first_node_
    buckets: Vec<usize>, // __bucket_list_: index of the node *before* the bucket's first node, or NIL
    size: usize,
    max_load_factor: f32,
}

impl<K: AsRef<[u8]>> Default for LibcxxHashSet<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: AsRef<[u8]>> LibcxxHashSet<K> {
    pub fn new() -> Self {
        LibcxxHashSet {
            nodes: vec![Node { next: NIL, hash: 0, value: None }],
            buckets: Vec::new(),
            size: 0,
            max_load_factor: 1.0,
        }
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// `bucket_count()`
    pub fn bucket_count(&self) -> usize {
        self.buckets.len()
    }

    fn find_node(&self, key: &[u8], hash: u64) -> Option<usize> {
        let bc = self.buckets.len();
        if bc == 0 {
            return None;
        }
        let chash = constrain_hash(hash, bc);
        let pn = self.buckets[chash];
        if pn == NIL {
            return None;
        }
        let mut nd = self.nodes[pn].next;
        while nd != NIL {
            let h = self.nodes[nd].hash;
            if !(h == hash || constrain_hash(h, bc) == chash) {
                break;
            }
            if h == hash && self.nodes[nd].value.as_ref().unwrap().as_ref() == key {
                return Some(nd);
            }
            nd = self.nodes[nd].next;
        }
        None
    }

    pub fn contains(&self, key: &[u8]) -> bool {
        self.find_node(key, libcxx_string_hash(key)).is_some()
    }

    /// `insert(value)`: returns true if newly inserted (`__emplace_unique_key_args`).
    pub fn insert(&mut self, value: K) -> bool {
        let hash = libcxx_string_hash(value.as_ref());
        if self.find_node(value.as_ref(), hash).is_some() {
            return false;
        }
        let mut bc = self.buckets.len();
        let nd = self.nodes.len();
        self.nodes.push(Node { next: NIL, hash, value: Some(value) });
        if ((self.size + 1) as f32) > (bc as f32) * self.max_load_factor {
            let a = 2 * bc + (!is_hash_power2(bc)) as usize;
            let b = (((self.size + 1) as f32) / self.max_load_factor).ceil() as usize;
            self.rehash_unique(a.max(b));
            bc = self.buckets.len();
        }
        let chash = constrain_hash(hash, bc);
        let pn = self.buckets[chash];
        if pn == NIL {
            let pn = BEFORE_BEGIN;
            self.nodes[nd].next = self.nodes[pn].next;
            self.nodes[pn].next = nd;
            self.buckets[chash] = pn;
            let nx = self.nodes[nd].next;
            if nx != NIL {
                let b = constrain_hash(self.nodes[nx].hash, bc);
                self.buckets[b] = nd;
            }
        } else {
            self.nodes[nd].next = self.nodes[pn].next;
            self.nodes[pn].next = nd;
        }
        self.size += 1;
        true
    }

    /// `__rehash<true>(n)`
    fn rehash_unique(&mut self, mut n: usize) {
        if n == 1 {
            n = 2;
        } else if n & (n - 1) != 0 {
            n = libcxx_next_prime(n as u64) as usize;
        }
        let bc = self.buckets.len();
        if n > bc {
            self.do_rehash(n);
        } else if n < bc {
            let need = ((self.size as f32) / self.max_load_factor).ceil() as usize;
            n = n.max(if is_hash_power2(bc) {
                next_hash_pow2(need)
            } else {
                libcxx_next_prime(need as u64) as usize
            });
            if n < bc {
                self.do_rehash(n);
            }
        }
    }

    /// `__do_rehash<true>(bucket_count)`
    fn do_rehash(&mut self, bucket_count: usize) {
        self.buckets.clear();
        self.buckets.resize(bucket_count, NIL);
        if bucket_count == 0 {
            return;
        }
        let mut pp = BEFORE_BEGIN;
        let mut cp = self.nodes[pp].next;
        if cp == NIL {
            return;
        }
        let mut chash = constrain_hash(self.nodes[cp].hash, bucket_count);
        self.buckets[chash] = pp;
        let mut phash = chash;
        pp = cp;
        cp = self.nodes[cp].next;
        while cp != NIL {
            chash = constrain_hash(self.nodes[cp].hash, bucket_count);
            if chash == phash {
                pp = cp;
            } else if self.buckets[chash] == NIL {
                self.buckets[chash] = pp;
                pp = cp;
                phash = chash;
            } else {
                let np = cp; // unique keys: no run of equal keys
                self.nodes[pp].next = self.nodes[np].next;
                let head = self.buckets[chash];
                self.nodes[np].next = self.nodes[head].next;
                self.nodes[head].next = cp;
            }
            cp = self.nodes[pp].next;
        }
    }

    /// Iterate in exactly the order of a C++ range-for over the set.
    pub fn iter(&self) -> LibcxxHashSetIter<'_, K> {
        LibcxxHashSetIter { set: self, cur: self.nodes[BEFORE_BEGIN].next }
    }
}

pub struct LibcxxHashSetIter<'a, K: AsRef<[u8]>> {
    set: &'a LibcxxHashSet<K>,
    cur: usize,
}

impl<'a, K: AsRef<[u8]>> Iterator for LibcxxHashSetIter<'a, K> {
    type Item = &'a K;
    fn next(&mut self) -> Option<&'a K> {
        if self.cur == NIL {
            return None;
        }
        let n = &self.set.nodes[self.cur];
        self.cur = n.next;
        n.value.as_ref()
    }
}

/// Emulates `std::unordered_set<std::string>` (default constructed, default
/// max_load_factor 1.0) populated only by `insert(std::string)`.
#[derive(Default)]
pub struct LibcxxStringSet {
    inner: LibcxxHashSet<String>,
}

impl LibcxxStringSet {
    pub fn new() -> Self {
        LibcxxStringSet { inner: LibcxxHashSet::new() }
    }
    /// Returns true if `s` was newly inserted.
    pub fn insert(&mut self, s: String) -> bool {
        self.inner.insert(s)
    }
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn bucket_count(&self) -> usize {
        self.inner.bucket_count()
    }
    pub fn contains(&self, s: &str) -> bool {
        self.inner.contains(s.as_bytes())
    }
    /// Elements in C++ iteration order.
    pub fn iter(&self) -> impl Iterator<Item = &str> + '_ {
        self.inner.iter().map(|s| s.as_str())
    }
    /// Elements in C++ iteration order, as bytes.
    pub fn iter_bytes(&self) -> impl Iterator<Item = &[u8]> + '_ {
        self.inner.iter().map(|s| s.as_bytes())
    }
}

// ---------------------------------------------------------------------------
// Boost 1.92 container_hash + delly hash helpers
// ---------------------------------------------------------------------------

#[inline]
fn boost_hash_mix(mut x: u64) -> u64 {
    const M: u64 = 0x0e9846af9b1a615d;
    x ^= x >> 32;
    x = x.wrapping_mul(M);
    x ^= x >> 32;
    x = x.wrapping_mul(M);
    x ^= x >> 28;
    x
}

#[inline]
fn mulx(x: u64, y: u64) -> u64 {
    let r = (x as u128) * (y as u128);
    (r as u64) ^ ((r >> 64) as u64)
}

#[inline]
fn read32le(s: &[u8], i: usize) -> u64 {
    u32::from_le_bytes(s[i..i + 4].try_into().unwrap()) as u64
}

/// `boost::hash_range(first, last)` for a random-access range of `char`
/// (e.g. `std::string_view` / `std::string` iterators, `const char*`), Boost 1.92, 64-bit.
/// Equivalent to `boost::hash<std::string>{}(s)` / `boost::hash_value(std::string)`.
pub fn boost_hash_range_bytes(s: &[u8]) -> u64 {
    boost_hash_range_bytes_seed(0, s)
}

/// `boost::hash_range(seed, first, last)` for a char range.
pub fn boost_hash_range_bytes_seed(seed: u64, s: &[u8]) -> u64 {
    const Q: u64 = 0x9e3779b97f4a7c15;
    const K: u64 = 0xdf442d22ce4859b9;
    let mut n = s.len();
    let mut p = 0usize;
    let mut w = mulx(seed.wrapping_add(Q), K);
    let mut h = w ^ (n as u64);
    while n >= 8 {
        let v1 = u64::from_le_bytes(s[p..p + 8].try_into().unwrap());
        w = w.wrapping_add(Q);
        h ^= mulx(v1.wrapping_add(w), K);
        p += 8;
        n -= 8;
    }
    let mut v1: u64 = 0;
    if n >= 4 {
        v1 = (read32le(s, p + n - 4) << ((n - 4) * 8)) | read32le(s, p);
    } else if n >= 1 {
        let x1 = (n - 1) & 2;
        let x2 = n >> 1;
        v1 = ((s[p + x1] as u64) << (x1 * 8)) | ((s[p + x2] as u64) << (x2 * 8)) | (s[p] as u64);
    }
    w = w.wrapping_add(Q);
    h ^= mulx(v1.wrapping_add(w), K);
    mulx(h.wrapping_add(w), K)
}

/// `boost::hash_combine(seed, v)` where `boost::hash_value(v)` is already computed.
#[inline]
pub fn boost_hash_combine_hashed(seed: &mut u64, hv: u64) {
    *seed = boost_hash_mix(seed.wrapping_add(0x9e3779b9).wrapping_add(hv));
}

/// `boost::hash_combine(seed, v)` with `v: std::size_t` / `uint64_t`
/// (`hash_value` is the identity).
pub fn boost_hash_combine_u64(seed: &mut u64, v: u64) {
    boost_hash_combine_hashed(seed, v);
}

/// `boost::hash_combine(seed, v)` with `v: int32_t` (`hash_value` = `static_cast<size_t>(v)`,
/// i.e. sign-extended).
pub fn boost_hash_combine_i32(seed: &mut u64, v: i32) {
    boost_hash_combine_hashed(seed, v as i64 as u64);
}

/// `boost::hash_combine(seed, v)` with `v: int64_t` (`hash_value` = `static_cast<size_t>(v)`).
pub fn boost_hash_combine_i64(seed: &mut u64, v: i64) {
    boost_hash_combine_hashed(seed, v as u64);
}

/// delly `hash_string(const char*)`. `s` must not contain NUL (C-string semantics: the
/// computation stops at the first NUL if present). `char` is signed on arm64 Apple, so
/// bytes >= 0x80 are sign-extended before the `* 76963` int multiply.
pub fn hash_string(s: &[u8]) -> u32 {
    let mut h: u32 = 37;
    for &c in s {
        if c == 0 {
            break;
        }
        let t = (c as i8 as i32).wrapping_mul(76963) as u32;
        h = h.wrapping_mul(54059) ^ t;
    }
    h
}

#[inline]
fn cstr(s: &[u8]) -> &[u8] {
    match s.iter().position(|&c| c == 0) {
        Some(i) => &s[..i],
        None => s,
    }
}

/// delly `hash_sr(bam1_t*)` given the read name and `core.flag`.
pub fn hash_sr(qname: &[u8], flag: u16) -> u64 {
    let q = cstr(qname);
    let mut seed = hash_string(q) as u64; // unsigned -> size_t: zero-extended
    let prehash = boost_hash_range_bytes(q);
    boost_hash_combine_u64(&mut seed, prehash);
    if (flag & 0x40) != 0 && seed > 0 {
        seed -= 1;
    }
    seed
}

/// delly `hash_pair(bam1_t*)`: combine order tid, pos, mtid, mpos.
pub fn hash_pair(qname: &[u8], tid: i32, pos: i64, mtid: i32, mpos: i64) -> u64 {
    let mut seed = hash_string(cstr(qname)) as u64;
    boost_hash_combine_i32(&mut seed, tid);
    boost_hash_combine_i64(&mut seed, pos);
    boost_hash_combine_i32(&mut seed, mtid);
    boost_hash_combine_i64(&mut seed, mpos);
    seed
}

/// delly `hash_pair_mate(bam1_t*)`: combine order mtid, mpos, tid, pos.
pub fn hash_pair_mate(qname: &[u8], tid: i32, pos: i64, mtid: i32, mpos: i64) -> u64 {
    let mut seed = hash_string(cstr(qname)) as u64;
    boost_hash_combine_i32(&mut seed, mtid);
    boost_hash_combine_i64(&mut seed, mpos);
    boost_hash_combine_i32(&mut seed, tid);
    boost_hash_combine_i64(&mut seed, pos);
    seed
}

#[cfg(test)]
mod tests {
    use super::*;

    // Golden values produced by C++ (Apple clang 21, libc++, -std=c++17 -O3 -DNDEBUG,
    // Boost 1.92) on arm64 macOS.
    const GOLDEN: &[(&[u8], u64, u64, u32)] = &[
        // (s, std::hash<std::string>, boost::hash_range, delly hash_string)
        (b"", 0x9ae16a3b2f90404f, 0x1c97da6f7f24eba6, 0x25),
        (b"a", 0x2420662cd003acfa, 0x5dc94fb3ff4ecd1c, 0x6f6cf4),
        (b"abc", 0x3a912f483a4ece31, 0x233d0900c00abcfd, 0x499f72d7),
        (b"ACGT", 0x7bf8a1619cadadc0, 0xf0fa2dd1b5570f37, 0xce081eba),
        (b"chrA_1_2", 0xa39fe0ffca1b8bb3, 0xb344806cd97c3a47, 0x85e7f5d6),
        (b"ACGTACGTAC", 0x7a8d954c488b47a8, 0xb6cdb58f52e13742, 0x8da0f2a1),
        (b"chrA_123_456_0:0:0_1:0:0_1a2b", 0xe9c6c33532d3c3d8, 0x0a97d3807232de28, 0x2c161262),
        (b"\xff\x80\x7f read", 0xb0996b35de1d89d3, 0x435b3157aa2f74ae, 0x0fb55511),
        (
            b"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ",
            0x53e75cf1c36b19bf,
            0x356e4f23e56908e7,
            0x0ae1efe2,
        ),
    ];

    #[test]
    fn string_hashes_match_cxx() {
        for &(s, h, b, hs) in GOLDEN {
            assert_eq!(libcxx_string_hash(s), h, "std::hash {:?}", s);
            assert_eq!(boost_hash_range_bytes(s), b, "boost::hash_range {:?}", s);
            assert_eq!(hash_string(s), hs, "hash_string {:?}", s);
        }
        let long: Vec<u8> = (0..100).map(|i| b'A' + ((i * 7) % 26) as u8).collect();
        assert_eq!(libcxx_string_hash(&long), 0x302b1c22bab97a47);
        assert_eq!(boost_hash_range_bytes(&long), 0x4b1ffc7c85408cee);
    }

    #[test]
    fn delly_hashes_match_cxx() {
        let q = b"chrA_123_456_0:0:0_1:0:0_1a2b";
        assert_eq!(hash_sr(q, 0), 0xfa14b2dcd0985df3);
        assert_eq!(hash_sr(q, 0x40), 0xfa14b2dcd0985df2);
        assert_eq!(hash_pair(q, 0, 12345, -1, -1), 0xe821cfa0e65c4eb0);
        assert_eq!(hash_pair_mate(q, 0, 12345, -1, -1), 0x12bd117a88426809);
        assert_eq!(hash_pair(q, 3, 249000000000, 21, -7), 0x628db375dce6b925);
    }

    #[test]
    fn next_prime_matches_libcxx() {
        assert_eq!(libcxx_next_prime(0), 0);
        assert_eq!(libcxx_next_prime(1), 2);
        assert_eq!(libcxx_next_prime(211), 211);
        assert_eq!(libcxx_next_prime(212), 223);
        assert_eq!(libcxx_next_prime(1595), 1597);
        assert_eq!(libcxx_next_prime(3195), 3203);
    }

    #[test]
    fn set_order_matches_cxx() {
        let ins = [
            "ACGT", "TTGA", "ACGT", "chr1", "x", "yy", "zzz",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "q", "w", "e", "r", "t", "y", "u", "i",
            "o", "p",
        ];
        let mut st = LibcxxStringSet::new();
        let added: Vec<bool> = ins.iter().map(|s| st.insert(s.to_string())).collect();
        assert!(!added[2]);
        assert_eq!(st.len(), 17);
        assert_eq!(st.bucket_count(), 23);
        let order: Vec<&str> = st.iter().collect();
        assert_eq!(
            order,
            [
                "o", "i", "u", "t", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "r", "yy", "e",
                "w", "q", "y", "zzz", "x", "chr1", "p", "TTGA", "ACGT",
            ]
        );
    }
}
