## An alternative covariance construction: Sigma = D * normalise(exp(A)) * D.
##
## An alternative to `sdcovsqrt2cov!`, selected at runtime by
## `ctsem_cov_expm!(true)` so both routes can be compared in one build with
## nothing else differing. The default is off, and with it off `sdcovsqrt2cov!`
## reproduces `constraincorsqrt1` to 1e-16.
##
## Why it exists: `constraincorsqrt1` is not onto the space of correlation
## matrices, and a gradient optimiser strands on it -- at k=12, 1 of 30 matched
## starts reached the known unique optimum of an unrestricted covariance problem,
## against 30 of 30 here. Measured on a fitted state-dependent DIFFUSION model
## (k=6, 60 subjects, 1500 rows, dev1): this route reaches a better likelihood,
## -6569.960321 against -6569.961397, identically from three starting vectors,
## and costs 1.36x the fit time.
##
## Layout is deliberately identical to the sd/correlation form, so no model spec
## changes: `mat`'s diagonal is the standard deviation and its lower triangle is
## the off-diagonal coordinate. What changes is the meaning of that coordinate --
## here it is the entry of a hollow symmetric `A`, and for a lone pair the
## resulting correlation is exactly `tanh(A[i,j])`, i.e. the coordinate is
## Fisher's z. Unlike the `constraincorsqrt1` coordinate it is unbounded, and the
## construction is permutation equivariant and onto.
##
## Scratch is cached per (element type, dimension, thread) rather than threaded
## through the filter's workspace, which would mean touching every construction
## site. The `::ExpmCovScratch{T,d}` assertion on the lookup is load-bearing and
## not decoration: without it every field read downstream is a boxed `Any`
## access, and that cost alone swamped the construction -- a real gradient came
## out at 2.85x rather than 1.18x, and 36x the primitive-level projection. The
## `Dict` lookup that remains is a few tens of nanoseconds against a construction
## of ~1 microsecond, so it is not worth removing.
##
## Like `_sdcovsqrt2cov_pullback!`, the reverse pass recomputes the forward
## quantities rather than storing them per row; the forward buffer is long
## overwritten by the time the reverse sweep arrives.

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
    g::Vector{T}
    yd::Vector{T}
    eb::ExpBuffer{T,D}
    fb::ExpFrechetBuffer{T,D}
    # Cache of exp(A), in the spirit of `DiscretizationCache`'s exp_table:
    # `_ekf_predict_step!` rebuilds the process noise on every row whether or not
    # DIFFUSION depends on the state, so for a constant DIFFUSION the same A
    # would otherwise be exponentiated once per row for no reason.
    #
    # It has to be a TABLE rather than one slot. The scratch is keyed by element
    # type and dimension, so T0VAR, DIFFUSION and MANIFESTVAR of the same size
    # all land here and a single slot thrashes -- measured: a constant-DIFFUSION
    # fit still paid ~890,000 exponentials with one slot, against ~1,590 for the
    # sd/correlation route.
    #
    # Keyed on each A's full contents rather than a hash, so an entry cannot go
    # stale: a hit costs up to `_EXPM_CACHE_SLOTS` comparisons of O(d^2) against
    # an O(d^3) Pade evaluation.
    Atab::Vector{Matrix{T}}
    Ytab::Vector{Matrix{T}}
    nslots::Base.RefValue{Int}
    nextslot::Base.RefValue{Int}
end

function ExpmCovScratch{T,D}() where {T,D}
    ExpmCovScratch{T,D}(zeros(T, D, D), zeros(T, D, D), zeros(T, D, D),
        zeros(T, D, D), zeros(T, D, D), zeros(T, D, D), zeros(T, D), zeros(T, D),
        ExpBuffer{T}(D), ExpFrechetBuffer{T}(D),
        [zeros(T, D, D) for _ in 1:_EXPM_CACHE_SLOTS],
        [zeros(T, D, D) for _ in 1:_EXPM_CACHE_SLOTS], Ref(0), Ref(1))
end

const _EXPM_COV_SCRATCH = Dict{Tuple{DataType,Int,Int},Any}()
const _EXPM_COV_LOCK = ReentrantLock()

# Keyed by thread as well as shape, so each thread works in its own buffers, and
# the first-use insertion is locked: without that, two threads reaching a new
# shape at once race on the Dict. The read is unlocked, which is safe here only
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

@inline function _blocks_equal(A::AbstractMatrix, B::AbstractMatrix, ::Val{d}) where {d}
    @inbounds for j in 1:d, i in 1:d
        A[i, j] == B[i, j] || return false
    end
    return true
end

# exp(A), with the content-keyed cache table in front of it
@inline function _cached_exp!(Y, sc, A, ::Val{d}) where {d}
    @inbounds for e in 1:sc.nslots[]
        if _blocks_equal(sc.Atab[e], A, Val(d))
            copyto!(Y, sc.Ytab[e])
            return nothing
        end
    end
    my_exp!(Y, A, sc.W1, sc.eb, Val(d))
    slot = sc.nextslot[]
    @inbounds copyto!(sc.Atab[slot], A)
    @inbounds copyto!(sc.Ytab[slot], Y)
    sc.nslots[] = max(sc.nslots[], slot)
    sc.nextslot[] = slot == _EXPM_CACHE_SLOTS ? 1 : slot + 1
    return nothing
end

# A is hollow symmetric, read from mat's lower triangle.
@inline function _fill_hollow!(A, mat, ::Val{d}, ::Type{T}) where {d,T}
    @inbounds for j in 1:d, i in 1:d
        A[i, j] = i == j ? zero(T) : (i > j ? T(mat[i, j]) : T(mat[j, i]))
    end
    return nothing
end

"""
    sdcovexpm2cov!(buffer, mat, dim)

Write `D * normalise(exp(A)) * D` into `buffer.out`, where `D` is `mat`'s
diagonal and `A` is hollow symmetric from `mat`'s lower triangle.
"""
function sdcovexpm2cov!(buffer, mat, ::Val{d}) where {d}
    T = eltype(buffer.out)
    sc = _expm_cov_scratch(T, Val(d))
    A, Y = sc.A, sc.Y
    _fill_hollow!(A, mat, Val(d), T)
    _cached_exp!(Y, sc, A, Val(d))
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
    A, Y = sc.A, sc.Y
    _fill_hollow!(A, mat, Val(d), T)
    _cached_exp!(Y, sc, A, Val(d))
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
    # Y = exp(A); the adjoint of E -> L(A, E) is Ybar -> L(A', Ybar), and A is
    # symmetric, so one Frechet evaluation at A in direction Ybar suffices.
    my_exp_frechet!(sc.Ytmp, sc.Abar, A, Ybar, sc.fb)
    @inbounds for i in 1:d, j in 1:(i - 1)
        mat_bar[i, j] += sc.Abar[i, j] + sc.Abar[j, i]
    end
    return mat_bar
end
