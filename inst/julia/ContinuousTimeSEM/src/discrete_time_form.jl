using ComponentArrays
using LinearAlgebra
using ForwardDiff

################################################################################
# Discrete-time EKF log-likelihood and gradient function
################################################################################

"""
Distinct `Δt` values an `ExpTable` holds for one `JAx`.

Sized for the shared wave schedule: a panel measured at the same irregular
times for every subject has as many distinct intervals as it has waves, and a
diary design with a handful of gap lengths has fewer still. Fully irregular
times overflow the table and replace entries round-robin, which costs a scan of
this many scalars per substep on top of the exponential they would have paid
anyway.
"""
const _CTSEM_EXP_TABLE = 32

"""
    ExpTable{T}(n)

The exponentials `exp(JAx * Δt)` seen for one `JAx`, keyed by `Δt`.

A panel design hands the filter the same `JAx` at every row of every subject
and a small set of intervals, in an order the last-value cache this replaces
could not exploit: a shared schedule of nineteen irregular waves missed on
every one of 3800 rows, because no row's interval equalled the previous one.
Here a lookup compares `JAx` once and then scans `dts[1:count]`.

The table is a separate mutable object from the `DiscretizationCache` that
holds it so that the subjects of one chunk can share a single table: the
reverse pass installs its chunk's table into each subject workspace before
filtering it (see `_ctsem_adjoint_chunk!`). One task runs a chunk, so that
sharing is race-free; two chunks never touch the same table.

A change of `JAx` -- a new trial point, a state-dependent model, a subject
whose drift carries a TI-predictor effect -- empties the table. Nothing here
assumes linearity; the `JAx` comparison is exact, and under `ForwardDiff.Dual`
it compares partials as well as values (see `_blocks_identical`).
"""
mutable struct ExpTable{T}
    JAx::Matrix{T}
    valid::Bool
    dts::Vector{T}
    outs::Vector{Matrix{T}}
    count::Int
    next::Int
end

function ExpTable(::Type{T}, n::Int, capacity::Int=_CTSEM_EXP_TABLE) where {T}
    return ExpTable{T}(zeros(T, n, n), false, zeros(T, capacity),
        [zeros(T, n, n) for _ in 1:capacity], 0, 1)
end

"""
    DiscretizationCache{T}

Cache for the two expensive, row-invariant pieces of
`_compute_discrete_time_form!`.

  * **`exp`**: `eJAx = exp(JAx * Δt)`, in an `ExpTable` keyed by `Δt` for the
    current `JAx`. It covers the balanced panel, the shared irregular schedule
    and the few-distinct-gaps diary design, and correctly misses on genuinely
    irregular observation times.
  * **`lyap`**: the asymptotic-diffusion solve `X` depends on `JAx` and the
    diffusion covariance and **not on `Δt` at all**, so for any model whose
    drift and diffusion are not state-dependent it is the same at every row of
    every subject. This is the important one: `schur!` was 38% of the primal
    filter's runtime and 60 of its 68 MB of allocations on a 20-latent model.
    A last-value entry is enough for it.

A state-dependent model changes `JAx` per row, so both simply miss every time
and pay only the comparison.
"""
mutable struct DiscretizationCache{T}
    exp_table::ExpTable{T}
    lyap_JAx::Matrix{T}
    lyap_Q::Matrix{T}
    lyap_out::Matrix{T}
    lyap_valid::Bool
end

function DiscretizationCache(::Type{T}, n::Int, k::Int) where {T}
    return DiscretizationCache{T}(ExpTable(T, n),
        zeros(T, k, k), zeros(T, k, k), zeros(T, k, k), false)
end

"""
    _blocks_identical(a, b, d1, d2)

Whether the leading `d1 x d2` blocks of `a` and `b` are identical.

Plain `==` is the right comparison here **including for `ForwardDiff.Dual`
entries**, which is worth stating because it is easy to assume otherwise:
ForwardDiff defines `==` on same-tag duals as
`==(value(x), value(y)) && ==(partials(x), partials(y))` (see `dual.jl`), so
two duals agreeing in value but not in derivative correctly compare unequal.
If that were value-only, these caches would reuse a matrix exponential whose
derivative belongs to a different point -- a silently wrong gradient on the
default ForwardDiff path. `test_adjoint_primitives.jl` pins that ForwardDiff
behaviour, so a change to it fails there rather than here.
"""
@inline function _blocks_identical(a, b, d1::Int, d2::Int)
    @inbounds for j in 1:d2, i in 1:d1
        a[i, j] == b[i, j] || return false
    end
    return true
end

@inline function _copy_block!(dest, src, d1::Int, d2::Int)
    @inbounds for j in 1:d2, i in 1:d1
        dest[i, j] = src[i, j]
    end
    return dest
end

# `nothing` disables caching, for the convenience constructor below and any
# caller that has no workspace to hang a cache on.
@inline _exp_cache_lookup(::Nothing, JAx, Δt, n::Int) = 0
@inline _exp_cache_store!(::Nothing, JAx, Δt, out, n::Int) = nothing
@inline _lyap_cache_hit(::Nothing, JAx, Q, k::Int) = false
@inline _lyap_cache_store!(::Nothing, JAx, Q, out, k::Int) = nothing

"""
The table slot holding `exp(JAx * Δt)`, or 0 on a miss.

`Δt` is compared as a scalar before anything else, so a schedule of irregular
intervals pays one comparison per live entry and never a matrix comparison
beyond the single `JAx` check.
"""
@inline function _exp_cache_lookup(cache::DiscretizationCache, JAx, Δt, n::Int)
    table = cache.exp_table
    if table.valid && _blocks_identical(table.JAx, JAx, n, n)
        @inbounds for e in 1:table.count
            if table.dts[e] == Δt
                _CTSEM_OPCOUNT.exp_cache_hit[] += 1
                return e
            end
        end
    end
    _CTSEM_OPCOUNT.exp_cache_miss[] += 1
    return 0
end

@inline function _exp_cache_store!(cache::DiscretizationCache, JAx, Δt, out, n::Int)
    table = cache.exp_table
    if !(table.valid && _blocks_identical(table.JAx, JAx, n, n))
        _copy_block!(table.JAx, JAx, n, n)
        table.valid = true
        table.count = 0
        table.next = 1
    end
    capacity = length(table.dts)
    if table.count < capacity
        e = table.count + 1
        table.count = e
    else
        e = table.next
        table.next = e == capacity ? 1 : e + 1
    end
    @inbounds table.dts[e] = Δt
    _copy_block!(table.outs[e], out, n, n)
    return nothing
end

@inline function _lyap_cache_hit(cache::DiscretizationCache, JAx, Q, k::Int)
    hit = cache.lyap_valid && _blocks_identical(cache.lyap_JAx, JAx, k, k) &&
        _blocks_identical(cache.lyap_Q, Q, k, k)
    (hit ? _CTSEM_OPCOUNT.lyap_cache_hit : _CTSEM_OPCOUNT.lyap_cache_miss)[] += 1
    return hit
end

@inline function _lyap_cache_store!(cache::DiscretizationCache, JAx, Q, out, k::Int)
    _copy_block!(cache.lyap_JAx, JAx, k, k)
    _copy_block!(cache.lyap_Q, Q, k, k)
    _copy_block!(cache.lyap_out, out, k, k)
    cache.lyap_valid = true
    return nothing
end

"""
    _compute_discrete_time_form!(discrete_ca, buffer, DIFFUSIONcov, pars, Δt, exp_buffer, lyap_buffer)

Compute continuous-to-discrete dynamics for one EKF time interval.

The function writes the discrete drift, intercept, and process covariance into
`discrete_ca` using the supplied scratch buffers. `DIFFUSIONcov` is the
continuous-time process covariance and `Δt` is the elapsed time.
"""
function _compute_discrete_time_form!(discrete_ca, buffer, DIFFUSIONcov, pars, Δt, exp_buffer, lyap_buffer)
    n = size(DIFFUSIONcov, 1)
    diffusion_state_indices = collect(1:n)
    diffusion_buffer = _make_square_buffer(eltype(DIFFUSIONcov), n)
    discretization_buffer = _make_discretization_buffer(eltype(DIFFUSIONcov), n)
    return _compute_discrete_time_form!(
        discrete_ca,
        buffer,
        DIFFUSIONcov,
        pars,
        Δt,
        exp_buffer,
        lyap_buffer,
        zeros(eltype(DIFFUSIONcov), n),
        diffusion_state_indices,
        diffusion_buffer,
        discretization_buffer,
        exp_buffer.dim,
        nothing,
    )
end

function _compute_discrete_time_form!(discrete_ca, buffer, DIFFUSIONcov, pars, Δt, exp_buffer,
    lyap_buffer, state, diffusion_state_indices, diffusion_buffer,
    discretization_buffer, dim::Val{d}, cache=nothing) where {d}
    # Stan propagates the local affine EKF model, not the raw DRIFT matrix:
    # f(x) = DRIFT * x + CINT, J = JAx, c = f(x) - J * x.
    slot = _exp_cache_lookup(cache, pars.JAx, Δt, d)
    if slot != 0
        _copy_block!(discrete_ca.eJAx, cache.exp_table.outs[slot], d, d)
    else
        copyto!(exp_buffer.As, pars.JAx)
        rmul!(exp_buffer.As, Δt)
        my_exp!(discrete_ca.eJAx, exp_buffer.As, buffer.intermediate, exp_buffer, dim)
        _exp_cache_store!(cache, pars.JAx, Δt, discrete_ca.eJAx, d)
    end
    copyto!(discrete_ca.dDRIFT, discrete_ca.eJAx)

    # The Lyapunov solve is restricted to the genuine diffusion block. Static
    # augmented coordinates have no diffusion and are excluded, as in Stan's
    # derrind subset; their covariance still propagates through full eJAx.
    fill!(discrete_ca.dDIFFUSION, zero(eltype(discrete_ca.dDIFFUSION)))
    kdim = length(diffusion_state_indices)
    fill!(discrete_ca.dINT, zero(eltype(discrete_ca.dINT)))

    # Stan forms the local affine correction only in the dynamic state block.
    # The full eJAx above still supplies the dynamic-from-static contribution;
    # keeping this solve in the dynamic block avoids inverting the singular
    # augmented Jacobian induced by static random-effect coordinates.
    @inbounds for i in 1:kdim
        ii = diffusion_state_indices[i]
        affine = pars.CINT[ii]
        for j in 1:d
            affine += (pars.DRIFT[ii, j] - pars.JAx[ii, j]) * state[j]
        end
        diffusion_buffer.r[i] = affine
    end
    @inbounds for j in 1:kdim, i in 1:kdim
        ii = diffusion_state_indices[i]
        jj = diffusion_state_indices[j]
        diffusion_buffer.intermediate[i, j] = pars.JAx[ii, jj]
    end
    @inbounds for i in 1:kdim
        ii = diffusion_state_indices[i]
        correction = -diffusion_buffer.r[i]
        for q in 1:kdim
            qq = diffusion_state_indices[q]
            correction += discrete_ca.eJAx[ii, qq] * diffusion_buffer.r[q]
        end
        diffusion_buffer.s[i] = correction
    end
    # `diffusion_buffer.dim` is `Val(kdim)` fixed at construction; building
    # `Val(kdim)` from the runtime length here was a dynamic dispatch and an
    # allocation on every prediction substep (Profile.Allocs, dev1).
    _solve_square_system!(diffusion_buffer.intermediate, diffusion_buffer.s,
        diffusion_buffer.piv, diffusion_buffer.dim)
    @inbounds for i in 1:kdim
        discrete_ca.dINT[diffusion_state_indices[i]] = diffusion_buffer.s[i]
    end

    @inbounds for j in 1:kdim, i in 1:kdim
        ii = diffusion_state_indices[i]
        jj = diffusion_state_indices[j]
        discretization_buffer.input[i, j] = pars.JAx[ii, jj]
        diffusion_buffer.intermediate[i, j] = DIFFUSIONcov[ii, jj]
    end
    dynamic_jacobian = view(discretization_buffer.input, 1:kdim, 1:kdim)
    # `diffusion_buffer.intermediate` holds the diffusion covariance here and
    # is overwritten immediately below, so the cache has to copy it now.
    if _lyap_cache_hit(cache, dynamic_jacobian, diffusion_buffer.intermediate, kdim)
        _copy_block!(diffusion_buffer.out, cache.lyap_out, kdim, kdim)
    else
        my_lyap!(diffusion_buffer.out, dynamic_jacobian, diffusion_buffer.intermediate, lyap_buffer)
        _lyap_cache_store!(cache, dynamic_jacobian, diffusion_buffer.intermediate,
            diffusion_buffer.out, kdim)
    end
    @inbounds for r in 1:kdim, i in 1:kdim
        value = zero(eltype(diffusion_buffer.intermediate))
        for q in 1:kdim
            value += discrete_ca.eJAx[diffusion_state_indices[i], diffusion_state_indices[q]] *
                diffusion_buffer.out[q, r]
        end
        diffusion_buffer.intermediate[i, r] = value
    end
    @inbounds for j in 1:kdim, i in 1:kdim
        value = diffusion_buffer.out[i, j]
        for r in 1:kdim
            value -= diffusion_buffer.intermediate[i, r] *
                discrete_ca.eJAx[diffusion_state_indices[j], diffusion_state_indices[r]]
        end
        discrete_ca.dDIFFUSION[diffusion_state_indices[i], diffusion_state_indices[j]] = value
    end
    return 
end

"""
    _compute_one_step_form!(discrete_ca, DIFFUSIONcov, pars, state,
                            diffusion_state_indices, dim)

The discrete-time case, where there is nothing to discretize.

A discrete-time model's DRIFT, CINT and DIFFUSION are already the one-step
quantities, so the three expensive pieces of the continuous form collapse: the
transition is `JAx` itself, the process noise is the diffusion covariance, and
the intercept is the local affine offset with no solve around it. `Δt` plays no
part -- a discrete model advances one step per row whatever the recorded
interval, which is also what Stan does.

The affine offset runs over every state, not only the diffusing ones: there
is no solve here to keep away from the singular augmented block, and a state
left out of it silently loses its CINT -- which is what an isolated
deterministic latent (no diffusion, no coupling to a diffusing state) falls
victim to, since `.ctJuliaDerrind()` excludes it.

That makes this function depend on DRIFT being the *true* one-step map on
every row, including the static coordinates a random effect augments the
state with. Their diagonal must be 1, not 0: a static state carries forward
whole. `.ctJuliaAugmentRandomEffects()` sets it, matching what
`ctJacobian()` does to the copy of DRIFT it builds JAx from, and the offset
then cancels to zero on those rows the way it should. When it was left at 0
the offset came out as `-x[i]` and cancelled the state the transition had
just carried, zeroing every random-effect coordinate at every step.
"""
function _compute_one_step_form!(discrete_ca, DIFFUSIONcov, pars, state,
    diffusion_state_indices, dim::Val{d}) where {d}
    copyto!(discrete_ca.eJAx, pars.JAx)
    copyto!(discrete_ca.dDRIFT, pars.JAx)

    @inbounds for i in 1:d
        affine = pars.CINT[i]
        for j in 1:d
            affine += (pars.DRIFT[i, j] - pars.JAx[i, j]) * state[j]
        end
        discrete_ca.dINT[i] = affine
    end

    fill!(discrete_ca.dDIFFUSION, zero(eltype(discrete_ca.dDIFFUSION)))
    kdim = length(diffusion_state_indices)
    @inbounds for j in 1:kdim, i in 1:kdim
        ii = diffusion_state_indices[i]
        jj = diffusion_state_indices[j]
        discrete_ca.dDIFFUSION[ii, jj] = DIFFUSIONcov[ii, jj]
    end
    return discrete_ca
end
