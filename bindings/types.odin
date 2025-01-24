package blis

Symmetric :: struct($T: typeid) {
    data: []T,
    m: int,
}

Matrix :: struct($T: typeid) {
    data: []T,
    m: int,
    n: int,
}

column :: #force_inline proc "contextless" (A: Matrix($E), j: int) -> []E {
    return A.data[j * A.m:(j + 1) * A.m]
}

QR :: struct($T: typeid) {
    Q: Matrix(T),
    R: Matrix(T),
}