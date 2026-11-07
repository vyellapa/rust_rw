// Exact port of libc++ `std::sort` (generic / non-branchless path) as shipped
// with Apple clang 21 (clang-2100.3.34.2), header:
//   /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1/__algorithm/sort.h
// plus the heap helpers it reaches through `__partial_sort`
// (__algorithm/partial_sort.h, make_heap.h, sort_heap.h, pop_heap.h,
// push_heap.h (`__sift_up`), sift_down.h (`__sift_down`, `__floyd_sift_down`)).
//
// Purpose: reproduce the *exact* final permutation libc++ produces for
// comparators that are strict weak orderings but not total orders (ties
// between distinguishable elements), so that a Rust port of a C++ tool can
// emit byte-identical output.
//
// ---------------------------------------------------------------------------
// Which libc++ code path this mirrors
// ---------------------------------------------------------------------------
// `std::sort(first, last[, comp])` -> `__sort_impl<_ClassicAlgPolicy>`:
//   1. `__debug_randomize_range` - a no-op unless the user defines
//      `_LIBCPP_DEBUG_RANDOMIZE_UNSPECIFIED_STABILITY` (nothing in the SDK
//      defines it; `-O3 -DNDEBUG` does not). Ignored here.
//   2. `__sort_dispatch`: the extern-template `__sort<__less<T>&, T*>` overloads
//      are only selected for built-in arithmetic `T` with `__less<>`, `less<T>`,
//      `less<>` or `ranges::less`. For a user struct the generic overload runs:
//      depth = 2 * __bit_log2(len), then
//      `__introsort<_, _Comp&, It, __use_branchless_sort<_Comp, It>>`.
//   3. `__use_branchless_sort` = contiguous iterator && cheap-to-copy &&
//      `is_arithmetic<value_type>` && comparator desugars to less/greater.
//      For a struct with `operator<` (called via `std::__less<void, void>`) or
//      any custom comparator on a struct, `is_arithmetic` is false, so it is
//      false: the classic introsort with `__partition_with_equals_on_right`
//      and the branchy `__sort3/4/5` is used. The bitset partition and the
//      `__cond_swap` sorting networks are never used for such types and are not
//      ported.
//   4. `__check_strict_weak_ordering_sorted` and the `__debug_less` wrapper
//      (`__comp_ref_type`) only exist in `_LIBCPP_HARDENING_MODE_DEBUG`. The
//      SDK's `__config_site` sets `_LIBCPP_HARDENING_MODE_DEFAULT 2`
//      (= `_LIBCPP_HARDENING_MODE_NONE`), and `-DNDEBUG` does not change it.
//      The fast/extensive modes only add assertions, which never move elements.
//      So the release algorithm is exactly what is ported below.
//
// Also, `__make_heap` has an `is_arithmetic<value_type>` special case
// (`__assume_both_children`). It is false for structs; this port implements
// the struct (false) variant only.
// CAVEAT: if this is used on a plain arithmetic `T` (e.g. `i32`) with a custom
// comparator, libc++ would take the `__assume_both_children` heap variant in
// the (rare) heapsort fallback, and with `std::less` the branchless path. The
// port matches libc++ only for non-arithmetic element types (structs), which
// is what it is meant for.
//
// ---------------------------------------------------------------------------
// Implementation notes
// ---------------------------------------------------------------------------
// libc++ moves elements out into temporaries (`value_type __t(iter_move(i))`)
// and shifts "holes". Here every such hole movement is done with
// `slice::swap`, keeping the temporary inside the slice at the current hole
// position. Each hole step `*hole = move(*src); hole = src;` becomes
// `swap(hole, src); hole = src;`, so the value that libc++ holds in a
// temporary always sits at the hole index, and comparisons against the
// temporary become comparisons against `v[hole]`. The comparator sees the same
// values in the same order, and the final arrangement is identical.
// For the partition routines the pivot temporary stays at `v[begin]`, which
// libc++ never reads or writes during the scan for a valid strict weak
// ordering. The final `*begin = move(*pivot_pos); *pivot_pos = move(pivot)`
// becomes `swap(begin, pivot_pos)`.
//
// Because only swaps are used, there is no `unsafe`, `T` needs neither `Copy`
// nor `Clone`, and a panicking comparator leaves the slice as a valid
// permutation of its input with nothing dropped twice or leaked.

/// Strict-weak "less than", the equivalent of a C++ `operator<`.
pub trait CppLess {
    fn lt(&self, other: &Self) -> bool;
}

/// `std::sort(v.begin(), v.end())` for a type with `operator<`
/// (libc++ calls it through `std::__less<void, void>`).
pub fn cpp_sort<T: CppLess>(v: &mut [T]) {
    cpp_sort_by(v, |a: &T, b: &T| a.lt(b));
}

/// `std::sort(v.begin(), v.end(), comp)` where `less(a, b)` == `comp(a, b)`.
/// `less` must be a strict weak ordering, like the C++ requirement.
pub fn cpp_sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut less: F) {
    // __sort_dispatch
    let len = v.len();
    if len == 0 {
        // "log(0) is undefined, so don't try computing the depth"
        return;
    }
    // __depth_limit = 2 * std::__bit_log2(__to_unsigned_like(__last - __first))
    let depth = 2 * (usize::BITS - 1 - len.leading_zeros()) as isize;
    introsort(v, 0, len, &mut less, depth, true);
}

// ------------------------------------------------------------------------
// Small sorting networks (non-branchless overloads)
// ------------------------------------------------------------------------

/// `__sort3` (the `!__use_branchless_sort` overload). Returns true if it swapped.
fn sort3<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x: usize, y: usize, z: usize, c: &mut F) -> bool {
    if !c(&v[y], &v[x]) {
        // x <= y
        if !c(&v[z], &v[y]) {
            return false;
        }
        v.swap(y, z);
        if c(&v[y], &v[x]) {
            v.swap(x, y);
        }
        return true;
    }
    if c(&v[z], &v[y]) {
        // x > y > z
        v.swap(x, z);
        return true;
    }
    v.swap(x, y);
    if c(&v[z], &v[y]) {
        v.swap(y, z);
    }
    true
}

/// `__sort4` (non-branchless).
fn sort4<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x1: usize, x2: usize, x3: usize, x4: usize, c: &mut F) {
    sort3(v, x1, x2, x3, c);
    if c(&v[x4], &v[x3]) {
        v.swap(x3, x4);
        if c(&v[x3], &v[x2]) {
            v.swap(x2, x3);
            if c(&v[x2], &v[x1]) {
                v.swap(x1, x2);
            }
        }
    }
}

/// `__sort5` (non-branchless).
fn sort5<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    x1: usize,
    x2: usize,
    x3: usize,
    x4: usize,
    x5: usize,
    c: &mut F,
) {
    sort4(v, x1, x2, x3, x4, c);
    if c(&v[x5], &v[x4]) {
        v.swap(x4, x5);
        if c(&v[x4], &v[x3]) {
            v.swap(x3, x4);
            if c(&v[x3], &v[x2]) {
                v.swap(x2, x3);
                if c(&v[x2], &v[x1]) {
                    v.swap(x1, x2);
                }
            }
        }
    }
}

// ------------------------------------------------------------------------
// Insertion sorts
// ------------------------------------------------------------------------

/// `__insertion_sort` (guarded) on `[first, last)`.
fn insertion_sort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    if first == last {
        return;
    }
    let mut i = first + 1;
    while i != last {
        let mut j = i - 1;
        if c(&v[i], &v[j]) {
            // value_type __t(iter_move(i)); the value of __t now lives at the hole `j`.
            let mut k = j;
            j = i;
            loop {
                // *__j = iter_move(__k); __j = __k;
                v.swap(j, k);
                j = k;
                // while (__j != __first && __comp(__t, *--__k))
                if j == first {
                    break;
                }
                k -= 1;
                if !c(&v[j], &v[k]) {
                    break;
                }
            }
            // *__j = std::move(__t);  (already at j)
        }
        i += 1;
    }
}

/// `__insertion_sort_unguarded` on `[first, last)`; requires `first > 0` and
/// `v[first - 1]` <= every element in the range.
fn insertion_sort_unguarded<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    if first == last {
        return;
    }
    let mut i = first + 1;
    while i != last {
        let mut j = i - 1;
        if c(&v[i], &v[j]) {
            let mut k = j;
            j = i;
            loop {
                v.swap(j, k);
                j = k;
                // while (__comp(__t, *--__k)) - no bounds check (sentinel at first-1).
                // An invalid comparator could run past index 0 here; in Rust that
                // panics on the index instead of reading out of bounds.
                k -= 1;
                if !c(&v[j], &v[k]) {
                    break;
                }
            }
        }
        i += 1;
    }
}

/// `__insertion_sort_incomplete`: sorts `[first, last)` but gives up after
/// `__limit = 8` element insertions, returning whether the range is sorted.
fn insertion_sort_incomplete<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    c: &mut F,
) -> bool {
    match last - first {
        0 | 1 => return true,
        2 => {
            if c(&v[last - 1], &v[first]) {
                v.swap(first, last - 1);
            }
            return true;
        }
        3 => {
            sort3(v, first, first + 1, last - 1, c);
            return true;
        }
        4 => {
            sort4(v, first, first + 1, first + 2, last - 1, c);
            return true;
        }
        5 => {
            sort5(v, first, first + 1, first + 2, first + 3, last - 1, c);
            return true;
        }
        _ => {}
    }
    let mut j = first + 2;
    sort3(v, first, first + 1, j, c);
    const LIMIT: u32 = 8;
    let mut count: u32 = 0;
    let mut i = j + 1;
    while i != last {
        if c(&v[i], &v[j]) {
            let mut k = j;
            j = i;
            loop {
                v.swap(j, k);
                j = k;
                if j == first {
                    break;
                }
                k -= 1;
                if !c(&v[j], &v[k]) {
                    break;
                }
            }
            count += 1;
            if count == LIMIT {
                // return ++__i == __last;
                return i + 1 == last;
            }
        }
        j = i;
        i += 1;
    }
    true
}

// ------------------------------------------------------------------------
// Partitioning (pivot is at v[begin] throughout; see module notes)
// ------------------------------------------------------------------------

/// `__partition_with_equals_on_right` on `[begin, end)`, pivot at `v[begin]`.
/// Returns (pivot position, already_partitioned).
fn partition_with_equals_on_right<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    c: &mut F,
) -> (usize, bool) {
    let mut first = begin;
    let mut last = end;
    // Find the first element greater or equal to the pivot (guarded by median-of-3).
    loop {
        first += 1;
        if !c(&v[first], &v[begin]) {
            break;
        }
    }
    // Find the last element less than the pivot.
    if begin == first - 1 {
        // while (__first < __last && !__comp(*--__last, __pivot));
        while first < last {
            last -= 1;
            if c(&v[last], &v[begin]) {
                break;
            }
        }
    } else {
        // Guarded.
        loop {
            last -= 1;
            if c(&v[last], &v[begin]) {
                break;
            }
        }
    }
    let already_partitioned = first >= last;
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if !c(&v[first], &v[begin]) {
                break;
            }
        }
        loop {
            last -= 1;
            if c(&v[last], &v[begin]) {
                break;
            }
        }
    }
    // Move the pivot to its correct position.
    let pivot_pos = first - 1;
    if begin != pivot_pos {
        v.swap(begin, pivot_pos);
    }
    (pivot_pos, already_partitioned)
}

/// `__partition_with_equals_on_left` on `[begin, end)`, pivot at `v[begin]`.
/// Returns the iterator one past the pivot (`__first`).
fn partition_with_equals_on_left<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    c: &mut F,
) -> usize {
    let mut first = begin;
    let mut last = end;
    if c(&v[begin], &v[last - 1]) {
        // Guarded.
        loop {
            first += 1;
            if c(&v[begin], &v[first]) {
                break;
            }
        }
    } else {
        // while (++__first < __last && !__comp(__pivot, *__first)) {}
        loop {
            first += 1;
            if first >= last || c(&v[begin], &v[first]) {
                break;
            }
        }
    }
    if first < last {
        // Always guarded thanks to median-of-three in __introsort.
        loop {
            last -= 1;
            if !c(&v[begin], &v[last]) {
                break;
            }
        }
    }
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if c(&v[begin], &v[first]) {
                break;
            }
        }
        loop {
            last -= 1;
            if !c(&v[begin], &v[last]) {
                break;
            }
        }
    }
    let pivot_pos = first - 1;
    if begin != pivot_pos {
        v.swap(begin, pivot_pos);
    }
    first
}

// ------------------------------------------------------------------------
// Introsort driver
// ------------------------------------------------------------------------

/// `__introsort<_AlgPolicy, _Compare, _RandomAccessIterator, /*_UseBitSetPartition=*/false>`
fn introsort<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    c: &mut F,
    mut depth: isize,
    mut leftmost: bool,
) {
    // Upper bound for using insertion sort for sorting.
    const LIMIT: usize = 24;
    // Lower bound for using Tuckey's ninther technique for median computation.
    const NINTHER_THRESHOLD: usize = 128;
    loop {
        let len = last - first;
        match len {
            0 | 1 => return,
            2 => {
                if c(&v[last - 1], &v[first]) {
                    v.swap(first, last - 1);
                }
                return;
            }
            3 => {
                sort3(v, first, first + 1, last - 1, c);
                return;
            }
            4 => {
                sort4(v, first, first + 1, first + 2, last - 1, c);
                return;
            }
            5 => {
                sort5(v, first, first + 1, first + 2, first + 3, last - 1, c);
                return;
            }
            _ => {}
        }
        // Use insertion sort if the length of the range is below the specified limit.
        if len < LIMIT {
            if leftmost {
                insertion_sort(v, first, last, c);
            } else {
                insertion_sort_unguarded(v, first, last, c);
            }
            return;
        }
        if depth == 0 {
            // Fallback to heap sort as Introsort suggests:
            // std::__partial_sort(__first, __last, __last, __comp)
            heap_sort_range(v, first, last, c);
            return;
        }
        depth -= 1;
        {
            let half_len = len / 2;
            if len > NINTHER_THRESHOLD {
                sort3(v, first, first + half_len, last - 1, c);
                sort3(v, first + 1, first + (half_len - 1), last - 2, c);
                sort3(v, first + 2, first + (half_len + 1), last - 3, c);
                sort3(v, first + (half_len - 1), first + half_len, first + (half_len + 1), c);
                v.swap(first, first + half_len);
            } else {
                sort3(v, first + half_len, first, last - 1, c);
            }
        }
        // If not leftmost and the pivot equals the element just left of the
        // range, everything equal to the pivot goes left and needs no sorting.
        if !leftmost && !c(&v[first - 1], &v[first]) {
            first = partition_with_equals_on_left(v, first, last, c);
            continue;
        }
        let (mut i, already_partitioned) = partition_with_equals_on_right(v, first, last, c);
        // If we were given a perfect partition, see if insertion sort is quick...
        if already_partitioned {
            let fs = insertion_sort_incomplete(v, first, i, c);
            if insertion_sort_incomplete(v, i + 1, last, c) {
                if fs {
                    return;
                }
                last = i;
                continue;
            } else if fs {
                i += 1;
                first = i;
                continue;
            }
        }
        // Sort the left partition recursively and the right partition with
        // tail recursion elimination.
        introsort(v, first, i, c, depth, leftmost);
        leftmost = false;
        i += 1;
        first = i;
    }
}

// ------------------------------------------------------------------------
// Heap sort fallback: __partial_sort(first, last, last) with middle == last
//   -> __partial_sort_impl: __make_heap(first, last); (empty loop over
//      [middle, last)); __sort_heap(first, last).
// ------------------------------------------------------------------------

/// `__partial_sort(__first, __last, __last, __comp)` on `[first, last)`.
fn heap_sort_range<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    if first == last {
        return;
    }
    make_heap(v, first, last, c);
    sort_heap(v, first, last, c);
}

/// `__make_heap` with `__assume_both_children == false` (non-arithmetic value_type).
fn make_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    let n = last - first;
    if n > 1 {
        // for (__start = (__n - 2) / 2; __start >= 0; --__start)
        let mut start = (n - 2) / 2;
        loop {
            sift_down(v, first, c, n, start);
            if start == 0 {
                break;
            }
            start -= 1;
        }
    }
}

/// `__sift_down<_AlgPolicy, /*__assume_both_children=*/false>`.
fn sift_down<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, c: &mut F, len: usize, mut start: usize) {
    let mut child = start;
    if len < 2 || (len - 2) / 2 < child {
        return;
    }
    child = 2 * child + 1;
    if child + 1 < len && c(&v[first + child], &v[first + child + 1]) {
        // right-child exists and is greater than left-child
        child += 1;
    }
    // check if we are in heap-order
    if c(&v[first + child], &v[first + start]) {
        return;
    }
    // value_type __top(iter_move(__first + __start)); __top travels with `start`.
    loop {
        v.swap(first + start, first + child);
        start = child;
        if (len - 2) / 2 < child {
            break;
        }
        child = 2 * child + 1;
        if child + 1 < len && c(&v[first + child], &v[first + child + 1]) {
            child += 1;
        }
        // } while (!__comp(__first[__child], __top));
        if c(&v[first + child], &v[first + start]) {
            break;
        }
    }
    // __first[__start] = std::move(__top);  (already there)
}

/// `__sort_heap`.
fn sort_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, mut last: usize, c: &mut F) {
    let mut n = last - first;
    while n > 1 {
        pop_heap(v, first, last, c, n);
        last -= 1;
        n -= 1;
    }
}

/// `__pop_heap`.
fn pop_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, mut last: usize, c: &mut F, len: usize) {
    if len > 1 {
        // __top = iter_move(__first); hole at __first. __top travels with the hole.
        let mut hole = floyd_sift_down(v, first, c, len);
        last -= 1;
        if hole == last {
            // *__hole = std::move(__top);  (already there)
        } else {
            // *__hole = iter_move(__last); *__last = std::move(__top);
            v.swap(hole, last);
            hole += 1;
            sift_up(v, first, hole, c, hole - first);
        }
    }
}

/// `__floyd_sift_down`: moves the hole from `first` down to a leaf, always
/// following the larger child. Returns the final hole index.
fn floyd_sift_down<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, c: &mut F, len: usize) -> usize {
    let mut hole = first;
    let mut child_i = first;
    let mut child: usize = 0;
    loop {
        child_i += child + 1;
        child = 2 * child + 1;
        if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
            // right-child exists and is greater than left-child
            child_i += 1;
            child += 1;
        }
        // swap __hole with its largest child
        v.swap(hole, child_i);
        hole = child_i;
        // if __hole is now a leaf, we're done
        if child > (len - 2) / 2 {
            return hole;
        }
    }
}

/// `__sift_up` (from push_heap.h) on `[first, last)` with `len = last - first`.
fn sift_up<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, mut last: usize, c: &mut F, mut len: usize) {
    if len > 1 {
        len = (len - 2) / 2;
        let mut ptr = first + len;
        last -= 1;
        if c(&v[ptr], &v[last]) {
            // __t = iter_move(__last); __t travels with `last`.
            loop {
                v.swap(last, ptr);
                last = ptr;
                if len == 0 {
                    break;
                }
                len = (len - 1) / 2;
                ptr = first + len;
                // } while (__comp(*__ptr, __t));
                if !c(&v[ptr], &v[last]) {
                    break;
                }
            }
        }
    }
}
