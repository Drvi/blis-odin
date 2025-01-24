package blis

import "core:math"

// DLARFG generates a real elementary reflector H of order n
@(require_results)
dlarfg :: #force_inline proc "c" (x: []f64) -> (tau: f64, alpha: f64) {
    n := len(x)
    alpha = x[0]
    if n <= 1 { return }

    xnorm: f64
    for i := 1; i < n; i += 1 {
        xnorm += x[i] * x[i]
    }
    if xnorm == 0 {
        return
    }

    xnorm = math.sqrt(xnorm)

    beta := -math.sign(alpha) * math.sqrt(alpha * alpha + xnorm * xnorm)
    safmin := f64(math.F64_MIN) / f64(math.F64_EPSILON)
    knt := 0
    if math.abs(beta) < safmin {
        rsafmn := f64(1) / safmin

        for math.abs(beta) < safmin {
            knt += 1
            // Compute x := x * (rsafmn)
            for i := 0; i < n; i += 1 {
                x[i] *= rsafmn
            }
            alpha *= rsafmn
            beta *= rsafmn
        }

        xnorm = 0
        for i := 1; i < n; i += 1 {
            xnorm += x[i] * x[i]
        }
        beta = -math.sign(alpha) * math.sqrt(alpha * alpha + xnorm * xnorm)
    }

    tau = (beta - alpha) / beta
    alpha = 1 / (alpha - beta)

    // Compute x := x * (1/alpha)
    for i := 0; i < n; i += 1 {
        x[i] *= alpha
    }

    alpha = beta
    x[0] = 1.0
    return
}

@(require_results)
slarfg :: #force_inline proc "c" (x: []f32) -> (tau: f32, alpha: f32) {
    n := len(x)
    alpha = x[0]
    if n <= 1 { return }

    xnorm: f32
    for i := 1; i < n; i += 1 {
        xnorm += x[i] * x[i]
    }
    if xnorm == 0 {
        return
    }

    xnorm = math.sqrt(xnorm)

    beta := -math.sign(alpha) * math.sqrt(alpha * alpha + xnorm * xnorm)
    safmin := f32(math.F32_MIN) / f32(math.F32_EPSILON)
    knt := 0
    if math.abs(beta) < safmin {
        rsafmn := f32(1) / safmin

        for math.abs(beta) < safmin {
            knt += 1
            // Compute x := x * (rsafmn)
            for i := 0; i < n; i += 1 {
                x[i] *= rsafmn
            }
            alpha *= rsafmn
            beta *= rsafmn
        }

        xnorm = 0
        for i := 1; i < n; i += 1 {
            xnorm += x[i] * x[i]
        }
        beta = -math.sign(alpha) * math.sqrt(alpha * alpha + xnorm * xnorm)
    }

    tau = (beta - alpha) / beta
    alpha = 1 / (alpha - beta)

    // Compute x := x * (1/alpha)
    for i := 0; i < n; i += 1 {
        x[i] *= alpha
    }

    alpha = beta
    x[0] = 1.0
    return
}
larfg :: proc{ slarfg, dlarfg }

// DLARF applies a real elementary reflector H to a real m by n matrix C
larf :: proc(v: []$E, tau: E, c: Matrix(E), work: []E) {
    if tau == 0 { return }
    // apply H to the left of C
    // w := C * v
    gemv(E(0.0), work, c, E(1.0), v, transa=trans.BLIS_TRANSPOSE)

    // C := C - tau * w * v^T
    ger(c, work, v, -tau)
    return
}

// We pass i in to subset the matrix and vector
_larf :: proc(v: []$E, tau: E, c: Matrix(E), work: []E, i: int) {
    if tau == 0 { return }
    // apply H to the left of C
    // w := C * v
    _work := work[i:]
    _gemv(E(0.0), _work, c, E(1.0), v, i, transa=trans.BLIS_TRANSPOSE)

    // C := C - tau * w * v^T
    ger(c, _work, v, -tau)
    return
}
_gemv :: #force_inline proc "contextless" (beta: $E, y: []E, A: Matrix(E), alpha: E, x: []E, i: int, transa: trans = trans.BLIS_NO_TRANSPOSE, conjx: conj = conj.BLIS_NO_CONJUGATE) {
    m := len(y)
    n := len(x)
    _alpha := alpha
    _beta := beta
    a_data = A.data[i * A.m + i:]
    // Note a_data is a slice of the original matrix and we're using A.m as the stride
    bli_gemv(transa, conjx, dim(m), dim(n), &_alpha, raw_data(a_data), inc(A.m), inc(1), raw_data(x), inc(1), &_beta, raw_data(y), inc(1))
    return
}


geqr2 :: proc(A: Matrix($E), tau: []E, work: []E) {
    m := A.m
    n := A.n
    if m < n { n = m }
    for i := 0; i < n; i += 1 {
        // Apply H to A[i:m, i]
        v := column(A, i)[i:]
        tau[i], alpha = larfg(v)
        _larf(v, tau[i], A, work, i)
    }
    return
}