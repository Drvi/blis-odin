package main

import blis "./bindings"
import "core:fmt"

main :: proc() {
    N := 8
    m := make([]f64, N*N)
    defer delete(m)

    // works with a slice
    fmt.println(m) // all zeroes
    alpha := 1.0
    blis.bli_dsetm(
        blis.conj.BLIS_NO_CONJUGATE,
        blis.doff(0),
        blis.diag.BLIS_NONUNIT_DIAG,
        blis.uplo.BLIS_DENSE,
        blis.dim(N),
        blis.dim(N),
        &alpha,
        raw_data(m),
        blis.inc(N),
        blis.inc(1),
    )
    fmt.println(m) // all ones
    fmt.println()

    // Wrap the slice that we just used
    t: blis.Symmetric(f64)
    t.data = m
    t.m = N

    fmt.println(t)
    // setsym :: proc(s: Symmetric($E), alpha: E) {
    //     alpha := alpha
    //     bli_dsetm(
    //         conj.BLIS_NO_CONJUGATE,
    //         doff(0),
    //         diag.BLIS_NONUNIT_DIAG,
    //         uplo.BLIS_DENSE,
    //         dim(s.m),
    //         dim(s.m),
    //         &alpha,
    //         raw_data(s.data),
    //         inc(s.m),
    //         inc(1),
    //     )
    //     return
    // }
    blis.setsym(t, 2.0) // segfaults
    fmt.println(t)
}
