package blis

import "core:c"
import "base:runtime"

when ODIN_OS == .Linux { foreign import blis "../blis/lib/zen3/libblis.a" }

Symmetric :: struct($T: typeid) {
    data: []T,
    m: int,
}

make_symmetric :: proc($T: typeid, m: int) -> (out: Symmetric(T), err: runtime.Allocator_Error) #optional_allocator_error {
    data := make([]T, m*m) or_return
    out.data = data
    out.m = m
    return
}

delete_symmetric :: proc(m: Symmetric($T)) {
    delete(m.data)
}

gint :: distinct c.int
dim :: distinct gint
inc :: distinct gint
doff :: distinct gint

trans :: enum c.int {
    BLIS_NO_TRANSPOSE      = 0,
    BLIS_TRANSPOSE         = 1 << 3,
    BLIS_CONJ_NO_TRANSPOSE =          1 << 4,
    BLIS_CONJ_TRANSPOSE    = 1 << 3 | 1 << 4,
}

conj :: enum c.int {
    BLIS_NO_CONJUGATE = 0,
    BLIS_CONJUGATE    = 1 << 4,
}

side :: enum c.int {
    BLIS_LEFT  = 0,
    BLIS_RIGHT = 1,
}

uplo :: enum c.int {
    BLIS_LOWER =          1 << 6 | 1 << 7,
    BLIS_UPPER = 1 << 5 | 1 << 6,
    BLIS_DENSE = 1 << 5 | 1 << 6 | 1 << 7,
}

diag :: enum c.int {
    BLIS_NONUNIT_DIAG = 0,
    BLIS_UNIT_DIAG    = 1 << 8,
}

@(default_calling_convention="c")
foreign blis {
    @(link_prefix="bli_")
    init :: proc() ---
    @(link_prefix="bli_")
    finalize :: proc() ---

    // level 0v
    bli_ssetv :: proc(conjalpha: conj, n: dim, alpha: ^f32, x: [^]f32, incx: inc) ---
    bli_dsetv :: proc(conjalpha: conj, n: dim, alpha: ^f64, x: [^]f64, incx: inc) ---

    // level 1v
    bli_saddv :: proc(conjx: conj, n: dim, x: [^]f32, incx: inc, y: [^]f32, incy: inc) ---
    bli_daddv :: proc(conjx: conj, n: dim, x: [^]f64, incx: inc, y: [^]f64, incy: inc) ---

    bli_sdotv :: proc(conjx: conj, conjy: conj, n: dim, x: [^]f32, incx: inc, y: [^]f32, incy: inc, rho: ^f32) ---
    bli_ddotv :: proc(conjx: conj, conjy: conj, n: dim, x: [^]f64, incx: inc, y: [^]f64, incy: inc, rho: ^f64) ---


    // level 1m
    bli_ssetm :: proc(conjalpha: conj, diagoffa: doff, diaga: diag, uploa: uplo, m: dim, n: dim, alpha: ^f32, a: [^]f32, rsa: inc, csa: inc) ---
    bli_dsetm :: proc(conjalpha: conj, diagoffa: doff, diaga: diag, uploa: uplo, m: dim, n: dim, alpha: ^f64, a: [^]f64, rsa: inc, csa: inc) ---

    //level 2
    bli_ssymv :: proc(uploa: uplo, conja: conj, conjx: conj, m: dim, alpha: ^f32, a: [^]f32, rsa: inc, csa: inc, x: [^]f32, incx: inc, beta: ^f32, y: [^]f32, incy: inc) ---
    bli_dsymv :: proc(uploa: uplo, conja: conj, conjx: conj, m: dim, alpha: ^f64, a: [^]f64, rsa: inc, csa: inc, x: [^]f64, incx: inc, beta: ^f64, y: [^]f64, incy: inc) ---
}

// SET
bli_setv :: proc {
    bli_ssetv,
    bli_dsetv,
}
bli_setm :: proc {
    bli_ssetm,
    bli_dsetm,
}

setv :: #force_inline proc "c" (x: $T/[]$E, alpha: E, n: int = -1, incx: int = 1, conjalpha: conj = conj.BLIS_NO_CONJUGATE) {
    _n := len(x) if n == -1 else n
    _alpha := alpha
    bli_setv(conjalpha, dim(_n), &_alpha, raw_data(x), inc(incx))
    return
}

setsym :: proc(s: Symmetric(f64), alpha: f64) {
    alpha := alpha
    bli_dsetm(
        conj.BLIS_NO_CONJUGATE,
        doff(0),
        diag.BLIS_NONUNIT_DIAG,
        uplo.BLIS_DENSE,
        dim(s.m),
        dim(s.m),
        &alpha,
        raw_data(s.data),
        inc(s.m),
        inc(1),
    )
    return
}

set :: proc {
    setv,
    setsym,
}

// ADD
bli_addv :: proc {
    bli_saddv,
    bli_daddv,
}

add :: #force_inline proc "c" (x, y: $T/[]$E, n: int = -1, incx: int = 1, incy: int = 1, conjx: conj = conj.BLIS_NO_CONJUGATE) {
    _n := len(x) if n == -1 else n
    bli_addv(conjx, dim(_n), raw_data(y), inc(incy), raw_data(x), inc(incx))
    return
}

// DOT
bli_dot :: proc {
    bli_sdotv,
    bli_ddotv,
}

@(require_results)
dot :: #force_inline proc "c" (x, y: $T/[]$E, n: int = -1, incx: int = 1, incy: int = 1, conjx: conj = conj.BLIS_NO_CONJUGATE, conjy: conj = conj.BLIS_NO_CONJUGATE) -> (rho: E) {
    _n := len(x) if n == -1 else n
    bli_dot(conjx, conjy, dim(_n), raw_data(x), inc(incx), raw_data(y), inc(incy), &rho)
    return
}


// SYMV
bli_symv :: proc {
    bli_ssymv,
    bli_dsymv,
}

symv :: #force_inline proc "c" (beta: $E, y: []E, alpha: E, a: Symmetric(E), x: []E, conja: conj = conj.BLIS_NO_CONJUGATE, conjx: conj = conj.BLIS_NO_CONJUGATE) {
    m := a.m
    _beta := beta
    _alpha := alpha

    bli_symv(uplo.BLIS_DENSE, conja, conjx, dim(m), &_alpha, raw_data(a.data), inc(m), inc(1), raw_data(x), inc(1), &_beta, raw_data(y), inc(1))
    return
}

