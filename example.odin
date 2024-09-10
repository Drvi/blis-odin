package main

import blis "./bindings"
import "core:fmt"


main :: proc() {
    blis.init()
    x := make([]f64, 10)
    defer delete(x)

    c := make([]f64, 10)
    defer delete(c)

    fmt.println(x)
    blis.set(x, 10.0)
    fmt.println(x)

    y := [?]f64{0.0, 0.0, 0.0}
    fmt.println(y)
    blis.set(y[:], 10.0)
    fmt.println(y)
    hadamard_product_v(y[:], y[:])
    fmt.println(y)
    z := [?]f64{1.0, 1.0, 1.0}
    blis.add(y[:], z[:])
    fmt.println(y)

    // fmt.println(blis.conj.BLIS_NO_CONJUGATE, " ", int(blis.conj.BLIS_NO_CONJUGATE), "\n",
    //             blis.diag.BLIS_NONUNIT_DIAG, " ", int(blis.diag.BLIS_NONUNIT_DIAG), "\n",
    //             blis.uplo.BLIS_DENSE, " ", int(blis.uplo.BLIS_DENSE))

    m, err := blis.make_symmetric(f64, 10)
    fmt.println(err, " ", len(m.data))
    // defer blis.delete_symmetric(m)
    fmt.println(m)
    // Segfaults without heap allocation
    d := new(f64)
    defer free(d)
    d^ = 2.0
    // d : ^f64
    // d^ = 2.0
    blis.set(m, d^)
    fmt.println(m)

    fmt.println(x)
    fmt.println(c)

    blis.symv(0.0, c, 1.0, m, x)
    fmt.println(x)
    fmt.println(c)

    rho := blis.dot(x, c)
    fmt.println(rho)

    defer blis.finalize()
}

hadamard_product_v :: #force_inline proc "c" (x, y: $T/[]$E) {
    n := len(x)
    for i in 0..<n {
        x[i] *= y[i]
    }
    return
}

// mutable struct OnlineRegressor
//     n::Int
//     coef::Vector{Float64}
//     Sxx_inv::Matrix{Float64}
//     Sxy::Vector{Float64}
//     OnlineRegressor(p) = new(0, zeros(p), zeros(p, p), zeros(p))
// end


// @inline function update_with_coefs!(o::OnlineRegressor, x::AbstractVector{T}, y::T, C::AbstractVector{T}) where T
//     o.n += 1
//     o.Sxy .+= y .* x

//     mul!(C, Symmetric(o.Sxx_inv), x)
//     s = 1 / (1 + x'C)
//     @inbounds @simd for i in axes(o.Sxx_inv, 1)
//         for j in axes(o.Sxx_inv, 2)
//             o.Sxx_inv[i, j] -= (C[i] * C[j]) * s
//         end
//     end

//     # if o.n > length(o.coef) # maybe slightly faster, how to initialize?
//     #     o.coef .= C .* ((y - x'o.coef) * s)
//     # else
//         mul!(o.coef, o.Sxx_inv, o.Sxy)
//     # end

//     return
// end
