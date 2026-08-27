"""
The joint posterior over population parameters *and* random effects, as a
log density with a gradient -- the target a Hamiltonian sampler needs.

# What this is an alternative to

`laplace.jl` integrates the random effects out by fitting a Gaussian at the mode
of each unit's inner problem. That is exact when the integrand is Gaussian in
`u` and approximate otherwise, and the error is not a constant: it grows
monotonically with the population scale, which tilts the profile and shrinks the
scale estimate. `quadrature.jl` measures that error and can correct it to first
order; this file removes it instead, by not making the approximation.

The target is

    log pi(theta, u) = sum_i loglik_i(theta + sum_l L_l(theta) u[b(i,l)])
                       - sum_U u_U'u_U / 2
                       + log p(theta)

which is the same object `_laplace_unit_objective_gradient` evaluates, summed
over units and with `theta` free rather than held. Sampling it gives the exact
posterior for both `theta` and every subject's effects, with no Gaussian
assumption anywhere.

# Why the parameterisation is already the right one

`u` enters as `theta + L(theta) u` with a standard normal prior, which is the
non-centred parameterisation -- the one hierarchical samplers need, because the
centred form couples each effect's scale to its location and produces Neal's
funnel. The Laplace layer chose it for its own reasons (the inner problem is
better conditioned that way) and the sampler inherits it for free.

It is not a complete escape: `L(theta)` still depends on the population scale,
so the funnel is weakened rather than removed, and a scale near zero still
pinches. `sample_nuts.jl` reports divergences, which is how that shows up.

# Cost

One gradient is one adjoint sweep per subject -- the same sweep the Laplace
inner Newton step uses, without the Newton iteration, without the curvature, and
without the log determinant. So a single gradient here is *cheaper* than a
single Laplace evaluation. The sampler needs far more of them, which is the
trade: exactness for time.

The dimension is `npar + sum_U dim(u_U)`, so a hundred subjects with two effects
each is a little over two hundred -- unremarkable for NUTS, and the reason the
metric in `sample_nuts.jl` is worth building carefully.
"""

using LinearAlgebra

export CTSEMSampler, ctsem_sampler, ctsem_sample_dimension

"""
    CTSEMSampler(laplace, npar)

The joint target built from a `CTSEMLaplaceObjective`.

Reuses the Laplace object outright -- its `spec` says which raw parameters vary
and at which level, its `units` say how the effects are laid out and shared, and
its per-chunk adjoint workspaces are exactly what the subject sweeps need. What
it does *not* reuse is `modes`: those are the Laplace approximation's answer,
and nothing here reads or writes them, so a fit may be sampled and then
summarised without the sample having moved anything underneath.

`x` is laid out as `[theta; u_1; u_2; ...]`, population parameters first and
then each unit's latent vector in unit order.
"""
struct CTSEMSampler{L}
    laplace::L
    npar::Int
    nunits::Int
    udims::Vector{Int}
    # Where unit `U`'s block of `x` begins, as a 0-based offset.
    uoffsets::Vector{Int}
    ndim::Int
    # Per level, the raw positions of its population scale and correlation
    # parameters -- the ones `L_l` depends on, and so the only ones with a
    # gradient contribution through the shift.
    positions::Vector{Vector{Int}}
end

function ctsem_sampler(laplace::CTSEMLaplaceObjective, npar::Integer)
    npar = Int(npar)
    npar >= 1 || throw(ArgumentError("npar must be positive"))
    _laplace_check_indices(laplace, npar)
    units = laplace.units
    nunits = length(units.members)
    udims = copy(units.dims)
    uoffsets = Vector{Int}(undef, nunits)
    at = npar
    for U in 1:nunits
        uoffsets[U] = at
        at += udims[U]
    end
    positions = [_laplace_level_positions(laplace.spec, l)
                 for l in eachindex(laplace.spec.levels)]
    return CTSEMSampler(laplace, npar, nunits, udims, uoffsets, at, positions)
end

"""Total dimension sampled: population parameters plus every unit's effects."""
ctsem_sample_dimension(sampler::CTSEMSampler) = sampler.ndim

"""Where unit `U`'s latent block sits in `x`, as a unit range."""
@inline _sample_urange(sampler::CTSEMSampler, U::Integer) =
    (sampler.uoffsets[U] + 1):(sampler.uoffsets[U] + sampler.udims[U])

"""
    ctsem_sample_density!(gradient, sampler, x)

The joint log density at `x`, writing its gradient into `gradient`.

Returns `-Inf` and leaves the gradient finite-but-meaningless when any subject's
likelihood is not finite, which a sampler must treat as a rejection rather than
an error: a leapfrog trajectory routinely steps somewhere the filter cannot
evaluate, and throwing there would end the chain instead of the trajectory.

Parallel over units, with the same cost-weighted chunking and the same
`ctsem_set_max_chunks!` ceiling the Laplace path uses. Each unit's latent block
is disjoint, so those gradients are written straight into `gradient`; only the
`theta` part needs a per-chunk accumulator, and it is `npar` long rather than
`ndim`.
"""
function ctsem_sample_density!(gradient::Vector{Float64}, sampler::CTSEMSampler,
    x::AbstractVector{Float64}; workspace_slot::Union{Nothing,Integer}=nothing)
    length(gradient) == sampler.ndim ||
        throw(DimensionMismatch("gradient must have $(sampler.ndim) entries"))
    length(x) == sampler.ndim ||
        throw(DimensionMismatch("x must have $(sampler.ndim) entries"))
    laplace = sampler.laplace
    spec = laplace.spec
    npar = sampler.npar
    nunits = sampler.nunits
    theta = Vector{Float64}(view(x, 1:npar))

    # The population Cholesky factors and their derivatives with respect to the
    # scale and correlation parameters. Both are `theta`-only, so they are built
    # once per gradient rather than once per subject.
    Ls, dL = try
        (_laplace_popchols(theta, spec), _laplace_level_chol_derivatives(theta, spec))
    catch err
        err isa InterruptException && rethrow()
        # A population covariance that will not factorize is outside the
        # support, not a bug: reject the point.
        fill!(gradient, 0.0)
        return -Inf
    end

    fill!(gradient, 0.0)
    # `workspace_slot` is how a *chain* claims a workspace. Several chains
    # running at once is the better parallel axis than several chunks within one
    # gradient -- chains share nothing and scale flat, where the unit loop
    # scales about twofold -- but they would then all be chunk 1 and race for
    # the same adjoint workspace. Given a slot, this runs serially inside it.
    parallel = workspace_slot === nothing
    nchunks = parallel ? _ctsem_nchunks(nunits) : 1
    slots = parallel ? (1:nchunks) : (Int(workspace_slot):Int(workspace_slot))
    while length(laplace.workspaces) < maximum(slots)
        push!(laplace.workspaces, Dict{Any,Any}())
    end
    ranges = parallel ?
        _ctsem_chunk_assignment(_laplace_unit_weights(laplace), nchunks) :
        [1:nunits]
    chunk_value = zeros(Float64, nchunks)
    chunk_theta = [zeros(Float64, npar) for _ in 1:nchunks]
    chunk_ok = fill(true, nchunks)

    run = function (c)
        slot = first(slots) + c - 1
        aws = _laplace_workspace!(laplace, Float64, npar, slot)
        # Per slot, not per subject: several chains filter the *same* subject at
        # the same time, where the unit loop never does, so the workspace cached
        # on the subject would be shared and silently corrupted.
        ekf = _laplace_ekf_workspace!(laplace, Float64, slot)
        gsub = Vector{Float64}(undef, npar)
        gtheta = chunk_theta[c]
        total = 0.0
        @inbounds for U in ranges[c]
            members = laplace.units.members[U]
            urange = _sample_urange(sampler, U)
            uview = view(x, urange)
            for (m, i) in enumerate(members)
                offsets = laplace.units.offsets[U][m]
                shifted = _laplace_member_values(theta, spec, Ls, uview, offsets)
                loglik = _laplace_subject_value_gradient!(gsub,
                    laplace.objective.subject_objectives[i], aws, shifted;
                    ekf_workspace=ekf)
                if !isfinite(loglik)
                    chunk_ok[c] = false
                    return nothing
                end
                total += loglik
                # `shifted = theta + sum_l L_l u`, so every parameter picks up
                # the subject's own gradient directly...
                for t in 1:npar
                    gtheta[t] += gsub[t]
                end
                for l in eachindex(spec.levels)
                    level = spec.levels[l]
                    k = nrandomeffects(level)
                    k == 0 && continue
                    base = offsets[l]
                    L = Ls[l]
                    re = level.re_index
                    # ...the effects pick it up through `L`...
                    for q in 1:k
                        acc = 0.0
                        for p in 1:k
                            acc += L[p, q] * gsub[re[p]]
                        end
                        gradient[sampler.uoffsets[U] + base + q] += acc
                    end
                    # ...and the scale and correlation parameters pick it up a
                    # second time, through `L` depending on them. Written as
                    # explicit loops rather than `dL[l][t] * u_block` because
                    # this is the innermost loop of the whole sampler and the
                    # product would allocate a vector per parameter per subject.
                    for (t, j) in enumerate(sampler.positions[l])
                        D = dL[l][t]
                        acc = 0.0
                        for p in 1:k
                            shift = 0.0
                            for q in 1:k
                                shift += D[p, q] * x[sampler.uoffsets[U] + base + q]
                            end
                            acc += gsub[re[p]] * shift
                        end
                        gtheta[j] += acc
                    end
                end
            end
            # The standard normal prior on this unit's effects.
            for a in urange
                total -= x[a] * x[a] / 2
                gradient[a] -= x[a]
            end
        end
        chunk_value[c] = total
        return nothing
    end

    if nchunks <= 1
        run(1)
    else
        Threads.@sync for c in 1:nchunks
            Threads.@spawn run(c)
        end
    end

    @inbounds for c in 1:nchunks
        chunk_ok[c] || return -Inf
    end
    value = sum(chunk_value)
    @inbounds for c in 1:nchunks, t in 1:npar
        gradient[t] += chunk_theta[c][t]
    end
    value += _ctsem_log_prior(laplace.objective, theta)
    _ctsem_log_prior_gradient!(view(gradient, 1:npar), laplace.objective, theta)
    return isfinite(value) ? value : -Inf
end

"""
    ctsem_sample_density(sampler, x)

Allocating form, for tests and for one-off checks.
"""
function ctsem_sample_density(sampler::CTSEMSampler, x::AbstractVector;
    workspace_slot::Union{Nothing,Integer}=nothing)
    g = zeros(Float64, sampler.ndim)
    value = ctsem_sample_density!(g, sampler, collect(Float64, x);
        workspace_slot=workspace_slot)
    return (value=value, gradient=g)
end

export ctsem_sample_density, ctsem_sample_density!

"""
    ctsem_sample_start(sampler, values)

An `x` holding the population parameters `values` and zero effects.

Zero is the prior mean of `u` and, after a Laplace fit, close to its posterior
mode as well -- the effects are standardised, so the modes are order one and the
chain is already in the typical set. That is the point of starting a sampler
from an optimised fit rather than from a random draw.
"""
function ctsem_sample_start(sampler::CTSEMSampler, values::AbstractVector)
    x = zeros(Float64, sampler.ndim)
    copyto!(view(x, 1:sampler.npar), view(collect(Float64, values), 1:sampler.npar))
    return x
end

export ctsem_sample_start
