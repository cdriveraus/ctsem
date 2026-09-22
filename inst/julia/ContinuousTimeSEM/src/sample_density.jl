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
its per-chunk adjoint workspaces are exactly what the subject sweeps need.

The joint density itself does not touch `modes`: the effects are coordinates of
`x` here, so a chain moves through them rather than solving for them, and a fit
can be sampled and then summarised without the sample having moved anything
underneath. Placing the chain's *starting* effects and metering their metric
blocks does read them, in `ctsem_sample_start` and `ctsem_sample_metric`, which
is why the latter solves them at the parameter vector it was given first.

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
    # Every subject of every unit, flattened, as parallel `(unit, member)`
    # vectors. The density loop divides over *this*, not over units: a sampled
    # model is frequently one study, and a unit axis then leaves every core but
    # one idle. Built once here because the alternative is rebuilding it on
    # every leapfrog step.
    flat_unit::Vector{Int}
    flat_member::Vector{Int}
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
    flat_unit = Int[]
    flat_member = Int[]
    for U in 1:nunits, m in eachindex(units.members[U])
        push!(flat_unit, U)
        push!(flat_member, m)
    end
    return CTSEMSampler(laplace, npar, nunits, udims, uoffsets, at, positions,
        flat_unit, flat_member)
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

Parallel over *subjects*, under the same `ctsem_set_max_chunks!` ceiling the
Laplace path uses.

Over units is what this did before, and it does not work: a sampled model is
often a single study, and there is then one unit and nothing to divide.
Measured on dev1 with 20 threads, one unit of 60 subjects, the same chain took
8.16 ms per leapfrog step on one core and 8.81 on twenty -- no speedup at all,
where the Laplace objective on the same shape gets 5.19x from dividing the same
subjects. The unit axis is the right one only when there are units to spare.

The price is that several subjects of one unit contribute to the *same* entries
of the gradient: the innermost level's block belongs to one subject, but every
level above it is shared by all the subjects under it, and `theta` is shared by
everyone. So each worker accumulates into a full-length buffer of its own and
adds it in once, under a lock, at the end of its chunk -- one acquisition per
worker per density call, against thousands of subject filters.

The standard normal prior on the effects is per unit, not per subject, so it
stays on the calling task above the loop.
"""
function ctsem_sample_density!(gradient::Vector{Float64}, sampler::CTSEMSampler,
    x::AbstractVector{Float64})
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
    # A chain is an item of the region one level out, so it already holds a
    # band of its own and splits it here; two chains can never reach the same
    # adjoint workspace because their bands are disjoint by construction.
    _laplace_ensure_pool!(laplace)

    # The standard normal prior on each unit's effects. Per unit, so it runs
    # here rather than inside the subject loop, and it touches each entry of
    # `gradient` exactly once.
    prior = 0.0
    @inbounds for U in 1:nunits
        for a in _sample_urange(sampler, U)
            prior -= x[a] * x[a] / 2
            gradient[a] -= x[a]
        end
    end

    nsub = length(sampler.flat_unit)
    reduction = ReentrantLock()
    shared_total = Ref(0.0)

    run = function (mine, _w)
        aws = _laplace_workspace!(laplace, Float64, npar)
        # Per slot, not per subject: several chains filter the *same* subject at
        # the same time, where the unit loop never does, so the workspace cached
        # on the subject would be shared and silently corrupted.
        ekf = _laplace_ekf_workspace!(laplace, Float64)
        gsub = _laplace_scratch_vector!(laplace, Float64, npar, :sample_gsub)
        shifted = _laplace_scratch_vector!(laplace, Float64, npar, :sample_shift)
        # This worker's whole contribution, added in once at the end. Full
        # length because subjects of one unit share every level above the
        # innermost, and all of them share `theta`.
        # `wacc`, not `acc`: the level loop below already has a scalar `acc`,
        # and a closure assigning a name the enclosing function also has is how
        # this engine loses a gradient silently.
        wacc = _laplace_scratch_vector!(laplace, Float64, sampler.ndim, :sample_acc)
        fill!(wacc, 0.0)
        gtheta = wacc
        total = 0.0
        @inbounds begin
            for s in mine
                U = sampler.flat_unit[s]
                m = sampler.flat_member[s]
                i = laplace.units.members[U][m]
                uview = view(x, _sample_urange(sampler, U))
                offsets = laplace.units.offsets[U][m]
                _laplace_member_values!(shifted, theta, spec, Ls, uview, offsets)
                loglik = _laplace_subject_value_gradient!(gsub,
                    laplace.objective.subject_objectives[i], aws, shifted;
                    ekf_workspace=ekf)
                if !isfinite(loglik)
                    return false
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
                        wacc[sampler.uoffsets[U] + base + q] += acc
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
        end
        Base.@lock reduction begin
            shared_total[] += total
            @inbounds for a in 1:sampler.ndim
                gradient[a] += wacc[a]
            end
        end
        return true
    end

    _laplace_partition(run, laplace, nsub) || return -Inf
    value = prior + shared_total[]
    value += _ctsem_log_prior(laplace.objective, theta)
    _ctsem_log_prior_gradient!(view(gradient, 1:npar), laplace.objective, theta)
    return isfinite(value) ? value : -Inf
end

"""
    ctsem_sample_density(sampler, x)

Allocating form, for tests and for one-off checks.
"""
function ctsem_sample_density(sampler::CTSEMSampler, x::AbstractVector;
    )
    g = zeros(Float64, sampler.ndim)
    value = ctsem_sample_density!(g, sampler, collect(Float64, x);
        )
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
function ctsem_sample_start(sampler::CTSEMSampler, values::AbstractVector;
    use_modes::Bool=true)
    # An entry point: the objective here is routinely a *different* one from
    # any this session has evaluated -- `ctSample()` on a reloaded fit builds a
    # fresh one -- and its workspace store is empty until this sizes it. The
    # band lives on the task and the store lives on the object, so arriving
    # with a band from some earlier object and no store of this one's is
    # exactly the mismatch `_laplace_check_slot` refuses.
    _laplace_ensure_pool!(sampler.laplace)
    x = zeros(Float64, sampler.ndim)
    copyto!(view(x, 1:sampler.npar), view(collect(Float64, values), 1:sampler.npar))
    use_modes || return x
    # The effects at their conditional modes rather than at zero.
    #
    # Zero is the prior mean and a defensible start, but the Laplace fit that
    # produced `values` has already solved for where each unit's effects
    # actually sit given those parameters, and that is a strictly better place
    # to begin: it is the mode of the very conditional distribution the metric's
    # effect blocks describe. Starting at zero asks the chain to travel there
    # first, through the part of warmup that is also adapting the step size.
    #
    # Guarded on length, and on nothing else -- in particular this cannot tell
    # a solved mode from an unsolved one. A `CTSEMLaplaceObjective` allocates
    # `modes` at full length and zero-filled, so an objective that has never
    # been evaluated reaches here with the guard satisfied and every effect at
    # zero. `ctsem_sample_metric` solves them at the parameter vector being
    # sampled from before this is called for exactly that reason; on a fit with
    # informative per-subject records, zero effects are not a neutral start but
    # a catastrophic one, measured at a joint density of -411899 where the fit's
    # own modes give -5297.
    laplace = sampler.laplace
    @inbounds for U in 1:sampler.nunits
        mode = laplace.modes[U]
        length(mode) == sampler.udims[U] || continue
        all(isfinite, mode) || continue
        for q in eachindex(mode)
            x[sampler.uoffsets[U] + q] = mode[q]
        end
    end
    return x
end

export ctsem_sample_start
