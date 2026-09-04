using LinearAlgebra

################################################################################
# Matrix exponential together with its Fréchet derivative
################################################################################

"""
    ExpFrechetBuffer{T}(n)

Workspace for `my_exp_frechet!`: the scaled inputs, the even powers of the
scaled matrix and their Fréchet counterparts, the Padé pieces for both, and two
scratch matrices. Eighteen `n × n` matrices in all.
"""
struct ExpFrechetBuffer{TYPE<:Number,D}
    As::Matrix{TYPE}
    Es::Matrix{TYPE}
    A2::Matrix{TYPE}
    A4::Matrix{TYPE}
    A6::Matrix{TYPE}
    M2::Matrix{TYPE}
    M4::Matrix{TYPE}
    M6::Matrix{TYPE}
    W1::Matrix{TYPE}
    W::Matrix{TYPE}
    Z1::Matrix{TYPE}
    V::Matrix{TYPE}
    U::Matrix{TYPE}
    LW::Matrix{TYPE}
    LU::Matrix{TYPE}
    LV::Matrix{TYPE}
    T1::Matrix{TYPE}
    T2::Matrix{TYPE}
    piv::Vector{Int}
    dim::Val{D}
    ExpFrechetBuffer{TYPE}(n::Int) where {TYPE<:Number} = new{TYPE,n}(
        (zeros(TYPE, n, n) for _ in 1:18)..., zeros(Int, n), Val(n))
end

# Al-Mohy & Higham (2009), Table 6.1: the largest norm at which the degree-13
# Padé approximant of the *Fréchet derivative* meets unit roundoff, which is a
# little tighter than the 5.37 `my_exp!` uses for the exponential alone.
const _FRECHET_THETA13 = 4.25

"""
    my_exp_frechet!(Y, L, A, E, buffer)

Write `exp(A)` to `Y` and the Fréchet derivative `L(A, E)` of the exponential at
`A` in direction `E` to `L`. `A` and `E` are read only.

This is the Al-Mohy & Higham (2009) recurrence: the Fréchet derivative of each
even power is carried alongside the power itself,

    M2 = A E + E A,   M4 = A2 M2 + M2 A2,   M6 = A4 M2 + M4 A2,

the Padé numerator and denominator are differentiated term by term, and the
quotient rule closes it: with `r = q⁻¹ p`, `L = q⁻¹ (Lp - Lq r)`. Squaring
then propagates both with `L ← R L + L R`, `R ← R²`. Everything is a product of
`n × n` matrices -- about nineteen of them plus two solves before squaring,
against six products of `2n × 2n` for the block identity `exp([A E; 0 A])`,
which is roughly eight times the arithmetic per product. It is exact, not an
approximation: both routes evaluate the same Padé approximant.

The direction enters linearly, so no normalisation of `E` is needed and none is
done; the scaling `s` is chosen from `‖A‖₁` alone, as it must be. A non-finite
input yields all-NaN output rather than an error, for the same reason `my_exp!`
does: an optimiser trial point that produces one is ordinary, and every caller
already treats a non-finite result as an invalid point.

The scalar type is generic, so `ForwardDiff.Dual` matrices go through the same
code with the engine's own LU -- which is what `ctsem_hessian` needs when it
differentiates the adjoint.
"""
function my_exp_frechet!(Y::AbstractMatrix{TYPE}, L::AbstractMatrix{TYPE},
    A::AbstractMatrix, E::AbstractMatrix, buf::ExpFrechetBuffer{TYPE,d}) where {TYPE<:Number,d}
    dim = buf.dim
    As, Es = buf.As, buf.Es
    @inbounds for j in 1:d, i in 1:d
        As[i, j] = A[i, j]
        Es[i, j] = E[i, j]
    end
    a1 = _opnorm1_noalloc(As, dim)
    e1 = _opnorm1_noalloc(Es, dim)
    if !(isfinite(a1) && isfinite(e1))
        fill!(Y, TYPE(NaN))
        fill!(L, TYPE(NaN))
        return Y, L
    end
    θ = _custom_abs(TYPE(_FRECHET_THETA13))
    s = a1 <= θ ? 0 : ceil(Int, log2(a1 / θ))
    factor = inv(TYPE(2)^s)
    rmul!(As, factor)
    rmul!(Es, factor)

    A2, A4, A6 = buf.A2, buf.A4, buf.A6
    M2, M4, M6 = buf.M2, buf.M4, buf.M6
    one_t = one(TYPE)
    mul!(A2, As, As)
    mul!(M2, As, Es)
    mul!(M2, Es, As, one_t, one_t)
    mul!(A4, A2, A2)
    mul!(M4, A2, M2)
    mul!(M4, M2, A2, one_t, one_t)
    mul!(A6, A4, A2)
    mul!(M6, A4, M2)
    mul!(M6, M4, A2, one_t, one_t)

    # b[1]..b[14] are b0..b13 of Higham (2008), as in `my_exp!`.
    b = (
        TYPE(64764752532480000.0), TYPE(32382376266240000.0), TYPE(7771770303897600.0),
        TYPE(1187353796428800.0), TYPE(129060195264000.0), TYPE(10559470521600.0),
        TYPE(670442572800.0), TYPE(33522128640.0), TYPE(1323241920.0), TYPE(40840800.0),
        TYPE(960960.0), TYPE(16380.0), TYPE(182.0), TYPE(1.0)
    )

    W1, W, Z1, V, U = buf.W1, buf.W, buf.Z1, buf.V, buf.U
    LW, LU, LV, T1, T2 = buf.LW, buf.LU, buf.LV, buf.T1, buf.T2

    # Odd part U = A (A6 W1 + W2), even part V = A6 Z1 + Z2.
    @. W1 = b[14] * A6 + b[12] * A4 + b[10] * A2
    @. Z1 = b[13] * A6 + b[11] * A4 + b[9] * A2
    mul!(W, A6, W1)
    @. W += b[8] * A6 + b[6] * A4 + b[4] * A2
    add_diag!(W, b[2], dim)
    mul!(V, A6, Z1)
    @. V += b[7] * A6 + b[5] * A4 + b[3] * A2
    add_diag!(V, b[1], dim)
    mul!(U, As, W)

    # Their Fréchet derivatives, term by term.
    @. T1 = b[14] * M6 + b[12] * M4 + b[10] * M2          # L(W1)
    mul!(LW, A6, T1)
    mul!(LW, M6, W1, one_t, one_t)
    @. LW += b[8] * M6 + b[6] * M4 + b[4] * M2            # L(W)
    mul!(LU, As, LW)
    mul!(LU, Es, W, one_t, one_t)                          # L(U) = A L(W) + E W
    @. T1 = b[13] * M6 + b[11] * M4 + b[9] * M2           # L(Z1)
    mul!(LV, A6, T1)
    mul!(LV, M6, Z1, one_t, one_t)
    @. LV += b[7] * M6 + b[5] * M4 + b[3] * M2            # L(V)

    # R = (V - U)⁻¹ (V + U), and L = (V - U)⁻¹ (LU + LV + (LU - LV) R).
    @. T1 = V - U
    @. Y = V + U
    copyto!(T2, T1)
    _solve_square_system!(T2, Y, buf.piv, dim)
    @. T2 = LU - LV
    mul!(L, T2, Y)
    @. L += LU + LV
    _solve_square_system!(T1, L, buf.piv, dim)

    # Undo the scaling: L(A², E) = A L + L A at each doubling.
    for _ in 1:s
        mul!(T1, Y, L)
        mul!(T1, L, Y, one_t, one_t)
        copyto!(L, T1)
        mul!(T1, Y, Y)
        copyto!(Y, T1)
    end
    return Y, L
end
