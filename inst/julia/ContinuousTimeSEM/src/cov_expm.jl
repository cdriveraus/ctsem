## An alternative covariance construction: Sigma = D * normalise(exp(A)) * D.
##
## Selected at runtime by `ctsem_cov_expm!(true)` so both routes can be compared
## in one build with nothing else differing. The default is off, and with it off
## `sdcovsqrt2cov!` reproduces `constraincorsqrt1` to 1e-16.
##
## Why it exists: `constraincorsqrt1` is not onto the space of correlation
## matrices, and a gradient optimiser strands on it -- at k=12, 1 of 30 matched
## starts reached the known unique optimum of an unrestricted covariance problem,
## against 30 of 30 here. It is also permutation equivariant, so an iid prior on
## the coordinates stays order invariant and a parameter shared across cells
## still means equal correlations. And its coefficients mean something: for a
## lone pair the correlation is exactly `tanh(A[i,j])`, so a coordinate is
## Fisher's z, unbounded rather than squashed into (-1, 1).
##
## Layout is deliberately identical to the sd/correlation form, so no model spec
## changes: `mat`'s diagonal is the standard deviation and its lower triangle is
## the off-diagonal coordinate, read into a hollow symmetric `A`.
##
## ## Cost, and where it goes
##
## `_ekf_predict_step!` calls `sdcovsqrt2cov!` unconditionally on every row, so
## the construction is rebuilt per row whether or not DIFFUSION depends on the
## state. Two things follow, and both are exploited below.
##
## `A` is hollow, so it carries only the correlation coordinates -- the standard
## deviations sit outside the exponential. In the usual formulation it is the
## variances that depend on the state, so `A` is then constant across rows and
## everything expensive can be cached against its contents.
##
## And `A` is symmetric, so the Frechet derivative needed by the reverse pass has
## a closed form through the eigendecomposition,
##
##     L(A, E) = Q [ (Q' E Q) .* Psi ] Q',
##     Psi[i,j] = (exp(l_i) - exp(l_j)) / (l_i - l_j),  exp(l_i) when l_i == l_j
##
## which is self-adjoint for symmetric `A`, so the same expression serves the
## pullback. `Q` and `Psi` depend on `A` alone, so a cache miss pays one
## symmetric eigendecomposition and each row then costs two matrix products
## instead of a full Pade recursion with its linear solve.
##
## The eigendecomposition needs a BLAS float, so anything else -- `Dual` above
## all -- falls back to `my_exp!` and `my_exp_frechet!`, which are correct at any
## element type. The fallback is exercised by the tests, not assumed.
##
## Scratch lives in a table keyed by (element type, dimension, thread). It has to
## be a table rather than one slot: T0VAR, DIFFUSION and MANIFESTVAR of the same
## size all share the entry, and with a single slot they evict each other -- a
## constant-DIFFUSION fit measured ~890,000 exponentials that way against the
## sd/correlation route's ~1,590. The `::ExpmCovScratch{T,d}` assertion on the
## lookup is load-bearing for the same reason: without it every field read is a
## boxed `Any` access, which alone put a real gradient at 2.85x rather than
## 1.18x.

const _CTSEM_COV_EXPM = Ref(false)

"""
    ctsem_cov_expm!(on::Bool)

Route `sdcovsqrt2cov!` and its pullback through the matrix-exponential
construction instead of `constraincorsqrt1`. Returns the previous setting.
"""
function ctsem_cov_expm!(on::Bool)
    prev = _CTSEM_COV_EXPM[]
    _CTSEM_COV_EXPM[] = on
    return prev
end

# Bang-free aliases: JuliaConnectoR addresses module members by name, and a
# trailing `!` is awkward to reach from R.
ctsem_cov_expm(on::Bool) = ctsem_cov_expm!(on)
ctsem_cov_expm() = _CTSEM_COV_EXPM[]

# Three distinct covariance matrices of one size is the common case (T0VAR,
# DIFFUSION, MANIFESTVAR); a few more slots cost only a comparison each.
const _EXPM_CACHE_SLOTS = 6

struct ExpmCovScratch{T,D}
    A::Matrix{T}
    Y::Matrix{T}
    W1::Matrix{T}
    Ybar::Matrix{T}
    Abar::Matrix{T}
    Ytmp::Matrix{T}
    M1::Matrix{T}
    M2::Matrix{T}
    g::Vector{T}
    yd::Vector{T}
    eb::ExpBuffer{T,D}
    fb::ExpFrechetBuffer{T,D}
    # Content-keyed cache. `Atab` holds the key, `Ytab` the exponential, and
    # `Qtab`/`Psitab` the eigendecomposition-derived Frechet factors when the
    # element type admits one (`eigok`).
    Atab::Vector{Matrix{T}}
    Ytab::Vector{Matrix{T}}
    Qtab::Vector{Matrix{T}}
    Psitab::Vector{Matrix{T}}
    eigok::Vector{Bool}
    nslots::Base.RefValue{Int}
    nextslot::Base.RefValue{Int}
end

function ExpmCovScratch{T,D}() where {T,D}
    mk() = [zeros(T, D, D) for _ in 1:_EXPM_CACHE_SLOTS]
    ExpmCovScratch{T,D}(zeros(T, D, D), zeros(T, D, D), zeros(T, D, D),
        zeros(T, D, D), zeros(T, D, D), zeros(T, D, D), zeros(T, D, D),
        zeros(T, D, D), zeros(T, D), zeros(T, D),
        ExpBuffer{T}(D), ExpFrechetBuffer{T}(D),
        mk(), mk(), mk(), mk(), fill(false, _EXPM_CACHE_SLOTS), Ref(0), Ref(1))
end

const _EXPM_COV_SCRATCH = Dict{Tuple{DataType,Int,Int},Any}()
const _EXPM_COV_LOCK = ReentrantLock()

# Keyed by thread as well as shape, so each thread works in its own buffers, and
# the first-use insertion is locked: without that, two threads reaching a new
# shape at once race on the Dict. The read is unlocked, which is safe only
# because an entry is never replaced once written -- a thread either sees
# `nothing` and takes the lock, or sees a fully constructed struct.
function _expm_cov_scratch(::Type{T}, ::Val{d}) where {T,d}
    key = (T, d, Threads.threadid())
    sc = get(_EXPM_COV_SCRATCH, key, nothing)
    if sc === nothing
        lock(_EXPM_COV_LOCK) do
            sc = get(_EXPM_COV_SCRATCH, key, nothing)
            if sc === nothing
                sc = ExpmCovScratch{T,d}()
                _EXPM_COV_SCRATCH[key] = sc
            end
        end
    end
    return sc::ExpmCovScratch{T,d}
end

# Only a BLAS float has a symmetric eigendecomposition here; `Dual` and friends
# take the Pade fallback.
_expm_eigen_eltype(::Type{T}) where {T<:Union{Float32,Float64}} = true
_expm_eigen_eltype(::Type) = false

@inline function _blocks_equal(A::AbstractMatrix, B::AbstractMatrix, ::Val{d}) where {d}
    @inbounds for j in 1:d, i in 1:d
        A[i, j] == B[i, j] || return false
    end
    return true
end

# A is hollow symmetric, read from mat's lower triangle.
@inline function _fill_hollow!(A, mat, ::Val{d}, ::Type{T}) where {d,T}
    @inbounds for j in 1:d, i in 1:d
        A[i, j] = i == j ? zero(T) : (i > j ? T(mat[i, j]) : T(mat[j, i]))
    end
    return nothing
end

"""
    _expm_cov_slot!(sc, A, dim)

Index of the cache slot holding `exp(A)` for this `A`, filling it on a miss.
"""
function _expm_cov_slot!(sc::ExpmCovScratch{T,D}, A, ::Val{d}) where {T,D,d}
    @inbounds for e in 1:sc.nslots[]
        _blocks_equal(sc.Atab[e], A, Val(d)) && return e
    end
    slot = sc.nextslot[]
    @inbounds begin
        copyto!(sc.Atab[slot], A)
        if _expm_eigen_eltype(T)
            # One symmetric eigendecomposition per distinct A. It allocates and
            # calls LAPACK, which is why it happens here and not per row.
            F = eigen(Symmetric(copy(A)))
            Q, lam = F.vectors, F.values
            copyto!(sc.Qtab[slot], Q)
            Psi = sc.Psitab[slot]
            for j in 1:d, i in 1:d
                delta = lam[i] - lam[j]
                # exp(l_j) * expm1(delta)/delta is the stable form of the
                # difference quotient, and tends to exp(l_j) as delta -> 0.
                Psi[i, j] = iszero(delta) ? exp(lam[j]) :
                    exp(lam[j]) * expm1(delta) / delta
            end
            # Y = Q diag(exp(lam)) Q', built through M1 to avoid a temporary.
            M1 = sc.M1
            for j in 1:d, i in 1:d
                M1[i, j] = Q[i, j] * exp(lam[j])
            end
            _mul_right_transpose!(sc.Ytab[slot], M1, Q, Val(d), Val(d), Val(d))
            sc.eigok[slot] = true
        else
            my_exp!(sc.Ytab[slot], A, sc.W1, sc.eb, Val(d))
            sc.eigok[slot] = false
        end
    end
    sc.nslots[] = max(sc.nslots[], slot)
    sc.nextslot[] = slot == _EXPM_CACHE_SLOTS ? 1 : slot + 1
    return slot
end

"""
    sdcovexpm2cov!(buffer, mat, dim)

Write `D * normalise(exp(A)) * D` into `buffer.out`, where `D` is `mat`'s
diagonal and `A` is hollow symmetric from `mat`'s lower triangle.
"""
function sdcovexpm2cov!(buffer, mat, ::Val{d}) where {d}
    T = eltype(buffer.out)
    sc = _expm_cov_scratch(T, Val(d))
    A = sc.A
    _fill_hollow!(A, mat, Val(d), T)
    slot = _expm_cov_slot!(sc, A, Val(d))
    Y = sc.Ytab[slot]
    yd, g = sc.yd, sc.g
    @inbounds for i in 1:d
        yd[i] = sqrt(Y[i, i])
        g[i] = T(mat[i, i]) / yd[i]
    end
    @inbounds for j in 1:d, i in 1:d
        buffer.out[i, j] = i == j ? T(mat[i, i])^2 : g[i] * g[j] * Y[i, j]
    end
    return nothing
end

# Frechet derivative of exp at A in direction E, through the cached
# eigendecomposition: L = Q [ (Q' E Q) .* Psi ] Q'.
@inline function _expm_frechet_eigen!(L, sc::ExpmCovScratch{T,D}, slot::Int, E,
        ::Val{d}) where {T,D,d}
    Q, Psi, M1, M2 = sc.Qtab[slot], sc.Psitab[slot], sc.M1, sc.M2
    # M2 = Q' E Q, elementwise-scaled by Psi. `mul!` reaches BLAS gemm here,
    # because this path only runs for a BLAS float.
    mul!(M1, transpose(Q), E)
    mul!(M2, M1, Q)
    @inbounds for j in 1:d, i in 1:d
        M2[i, j] *= Psi[i, j]
    end
    # L = Q * M2 * Q'
    mul!(M1, Q, M2)
    _mul_right_transpose!(L, M1, Q, Val(d), Val(d), Val(d))
    return nothing
end

"""
    _sdcovexpm2cov_pullback!(mat_bar, mat, cov_bar, d)

Reverse pass of `sdcovexpm2cov!`, accumulating into `mat_bar` in `mat`'s own
layout. Follows `_sdcovsqrt2cov_pullback!`'s convention: the cotangent is
symmetrised by adding both triangles.
"""
function _sdcovexpm2cov_pullback!(mat_bar::AbstractMatrix, mat::AbstractMatrix,
        cov_bar::AbstractMatrix, d::Int)
    d == 0 && return mat_bar
    all(iszero, cov_bar) && return mat_bar
    # One dynamic dispatch per call to reach the specialised body, rather than a
    # boxed field access for every read inside it.
    return _sdcovexpm2cov_pullback!(mat_bar, mat, cov_bar, Val(d))
end

function _sdcovexpm2cov_pullback!(mat_bar::AbstractMatrix, mat::AbstractMatrix,
        cov_bar::AbstractMatrix, ::Val{d}) where {d}
    T = promote_type(eltype(mat), eltype(cov_bar))
    sc = _expm_cov_scratch(T, Val(d))
    A = sc.A
    _fill_hollow!(A, mat, Val(d), T)
    slot = _expm_cov_slot!(sc, A, Val(d))
    Y = sc.Ytab[slot]
    yd = sc.yd
    @inbounds for i in 1:d
        yd[i] = sqrt(Y[i, i])
    end

    Ybar = sc.Ybar
    fill!(Ybar, zero(T))
    # Sigma_ii = s_i^2 exactly, so the diagonal carries no dependence on Y.
    # Off the diagonal Sigma_ij = s_i s_j C_ij with C_ij = Y_ij / (y_i y_j).
    @inbounds for i in 1:d
        acc = 2 * T(cov_bar[i, i]) * T(mat[i, i])
        for j in 1:d
            j == i && continue
            w = T(cov_bar[i, j]) + T(cov_bar[j, i])
            # For a fixed i this visits each j once, and `w` already carries
            # both of that pair's cotangent entries, so there is no halving.
            acc += w * T(mat[j, j]) * (Y[i, j] / (yd[i] * yd[j]))
        end
        mat_bar[i, i] += acc
    end
    # C -> Y. Only the symmetrisation of Ybar matters downstream, because the
    # Frechet map at a symmetric A sends symmetric directions to symmetric
    # results, so putting the whole pair weight in the lower triangle is exact.
    @inbounds for i in 1:d, j in 1:(i - 1)
        w = T(cov_bar[i, j]) + T(cov_bar[j, i])
        cbar = w * T(mat[i, i]) * T(mat[j, j])
        cij = Y[i, j] / (yd[i] * yd[j])
        Ybar[i, j] += cbar / (yd[i] * yd[j])
        Ybar[i, i] -= cbar * cij / (2 * Y[i, i])
        Ybar[j, j] -= cbar * cij / (2 * Y[j, j])
    end
    # Y = exp(A), and the adjoint of E -> L(A, E) is L(A', .) = L(A, .) for
    # symmetric A, so the same map serves the pullback.
    if sc.eigok[slot]
        _expm_frechet_eigen!(sc.Abar, sc, slot, Ybar, Val(d))
    else
        my_exp_frechet!(sc.Ytmp, sc.Abar, A, Ybar, sc.fb)
    end
    @inbounds for i in 1:d, j in 1:(i - 1)
        mat_bar[i, j] += sc.Abar[i, j] + sc.Abar[j, i]
    end
    return mat_bar
end
