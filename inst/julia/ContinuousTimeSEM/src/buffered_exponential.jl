using LinearAlgebra

################################################################################
# Matrix exponential related functions
################################################################################
"""
    ExpBuffer{T}(n)

Workspace for allocation-conscious matrix exponential computations.

The buffer stores scaled powers of the input matrix and pivot storage reused by
`my_exp!`.
"""
struct ExpBuffer{TYPE<:Number,D}
    As::Matrix{TYPE}
    A2::Matrix{TYPE}
    A4::Matrix{TYPE}
    A6::Matrix{TYPE}
    piv::Vector{Int}
    dim::Val{D}
    ExpBuffer{TYPE}(n::Int) where {TYPE <: Number} = new{TYPE,n}(
        zeros(TYPE, n, n), zeros(TYPE, n, n), zeros(TYPE, n, n), zeros(TYPE, n, n),
        zeros(Int, n),
        Val(n),
    )

end

"""
    ExpBuffer(mat::AbstractMatrix)

Create an `ExpBuffer` sized and typed from `mat`.
"""
ExpBuffer(mat::AbstractMatrix) = ExpBuffer{eltype(mat)}(size(mat, 1))

function _opnorm1_noalloc(A::AbstractMatrix{T}, ::Val{d}) where {T<:Number, d}
    z = zero(_custom_abs(zero(T)))
    max_col_sum = z
    @inbounds for j in 1:d
        col_sum = z
        for i in 1:d
            col_sum += _custom_abs(A[i, j])
        end
        max_col_sum = max(max_col_sum, col_sum)
    end
    return max_col_sum
end


"""
    my_exp!(Y, A, W1, buffer)

Compute the matrix exponential `exp(A)` and write it to `Y`.

The implementation uses Higham's scaling-and-squaring method with a `[13/13]`
Pade approximant. `W1` and `buffer` are scratch storage and may be overwritten.
"""
function my_exp!(Y::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, W1::AbstractMatrix{TYPE}, buffer::ExpBuffer{TYPE}) where {TYPE<:Number}
    return my_exp!(Y, A, W1, buffer, buffer.dim)
end

function my_exp!(Y::AbstractMatrix{TYPE}, A::AbstractMatrix{TYPE}, W1::AbstractMatrix{TYPE}, buffer::ExpBuffer{TYPE}, dim::Val{d}) where {TYPE<:Number, d}
    # Higham (2008) scaling-and-squaring with [13/13] Pade approximant.
    a1 = _opnorm1_noalloc(A, dim)
    # A matrix with a non-finite entry has no exponential, and saying so with
    # NaN is the only way to say it that the callers can act on.
    #
    # Without this the squaring count `ceil(Int, log2(a1/θ13))` is `Int(NaN)`,
    # which throws an `InexactError` from six frames inside the reverse pass.
    # An optimizer trial point that produces a non-finite DRIFT is an ordinary
    # event -- it is what a line search is for -- and every caller already
    # treats a non-finite objective as an invalid point and shrinks the step.
    # A thrown error instead ends the whole fit, and it ended it in R, where
    # `InexactError: Int64(NaN)` names neither the matrix nor the parameter
    # that produced it. Seen on one starting draw in ten of a 25-subject
    # Laplace fit.
    if !isfinite(a1)
        fill!(Y, TYPE(NaN))
        return Y
    end
    θ13 = _custom_abs(TYPE(5.371920351148152))
    s = a1 <= θ13 ? 0 : ceil(Int, log2(a1 / θ13))
    factor = inv(TYPE(2)^s)

    # Compute powers of A
    @inbounds for j in 1:d, i in 1:d
        buffer.As[i, j] = factor * A[i, j]
    end
    mul!(buffer.A2, buffer.As, buffer.As)
    mul!(buffer.A4, buffer.A2, buffer.A2)
    mul!(buffer.A6, buffer.A4, buffer.A2)

    # Coefficients b[1]..b[14] correspond to b0..b13 in Higham (2008).
    b = (
        TYPE(64764752532480000.0), TYPE(32382376266240000.0), TYPE(7771770303897600.0),
        TYPE(1187353796428800.0), TYPE(129060195264000.0), TYPE(10559470521600.0),
        TYPE(670442572800.0), TYPE(33522128640.0), TYPE(1323241920.0), TYPE(40840800.0),
        TYPE(960960.0), TYPE(16380.0), TYPE(182.0), TYPE(1.0)
    )

    # Odd polynomial core (before left-multiplication by As):
    # W = A6*(b13*A6 + b11*A4 + b9*A2) + b7*A6 + b5*A4 + b3*A2 + b1*I
    @. W1 = (b[14] * buffer.A6 + b[12] * buffer.A4 + b[10] * buffer.A2)
    mul!(Y, buffer.A6, W1)
    @. Y += b[8] * buffer.A6 + b[6] * buffer.A4 + b[4] * buffer.A2
    add_diag!(Y, b[2], dim)

    # U = As * W. Store U in W1.
    mul!(W1, buffer.As, Y)

    # Even polynomial:
    # V = A6*(b12*A6 + b10*A4 + b8*A2) + b6*A6 + b4*A4 + b2*A2 + b0*I
    @. Y = (b[13] * buffer.A6 + b[11] * buffer.A4 + b[9] * buffer.A2)
    mul!(buffer.As, buffer.A6, Y)
    @. buffer.As += b[7] * buffer.A6 + b[5] * buffer.A4 + b[3] * buffer.A2
    add_diag!(buffer.As, b[1], dim)

    # Form P = V + U (in Y) and Q = V - U (in A2), then solve Q * Y = P.
    @inbounds for j in 1:d, i in 1:d
        vij = buffer.As[i, j]
        uij = W1[i, j]
        buffer.A2[i, j] = vij - uij
        Y[i, j] = vij + uij
    end

    # Solve Q * Y = P in-place using reusable pivot storage.
    _solve_square_system!(buffer.A2, Y, buffer.piv, dim)

    # Squaring phase
    for _ in 1:div(s, 2)
    # for _ in 1:s
        mul!(W1, Y, Y)
        mul!(Y, W1, W1)
    end
    if isodd(s)
        mul!(W1, Y, Y)
        copyto!(Y, W1)
    end
    return Y
end
