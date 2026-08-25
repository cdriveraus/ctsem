using DataFrames, ForwardDiff, LinearAlgebra

# Laplace-approximate marginal likelihood for subject-level random effects.
#
# Two things need proving, and they need proving separately.
#
#  1. That the *approximation* is assembled correctly -- right mode, right
#     curvature, right constants. This is checked where the approximation is
#     not an approximation at all: when a random effect enters the state mean
#     linearly through an identity transform in a linear-Gaussian model, the
#     integrand is exactly Gaussian in `z`, Laplace is exact, and the integral
#     has a closed form. `_gaussian_reference` computes that closed form from
#     the log likelihood itself, by recovering the quadratic it must be, so the
#     reference shares no code with the thing it is checking.
#
#  2. That the *outer gradient* is the gradient of the value actually returned.
#     That is finite differences of the Laplace value, and it has to hold for
#     nonlinear models too, where the value is genuinely approximate. Getting
#     the approximation slightly wrong but differentiating it consistently is a
#     failure mode that the exactness check above cannot see, and vice versa.

################################################################################
# Model construction
################################################################################

function _laplace_test_dataframe(;
    drift::AbstractMatrix, jax::AbstractMatrix, cint::AbstractMatrix,
    diffusion::AbstractMatrix, lambda::AbstractMatrix, jy::AbstractMatrix,
    manifestmeans::AbstractMatrix, manifestvar::AbstractMatrix,
    t0var::AbstractMatrix, t0means::AbstractMatrix,
    pars::AbstractMatrix=zeros(1, 1),
    free::Dict=Dict(), predict::Dict=Dict(),
)
    matrices = Symbol[]; rows = Int[]; cols = Int[]
    parnumber = Union{Missing,Int}[]; value = Union{Missing,Float64}[]
    transform = Union{Missing,String}[]; predicttransform = Union{Missing,String}[]

    function addmat!(name::Symbol, mat::AbstractMatrix)
        for j in axes(mat, 2), i in axes(mat, 1)
            key = (name, i, j)
            push!(matrices, name); push!(rows, i); push!(cols, j)
            if haskey(free, key)
                pn, tf = free[key]
                push!(parnumber, pn); push!(value, missing); push!(transform, tf)
            else
                push!(parnumber, missing); push!(value, Float64(mat[i, j])); push!(transform, missing)
            end
            push!(predicttransform, get(predict, key, missing))
        end
    end

    addmat!(:DRIFT, drift); addmat!(:JAx, jax); addmat!(:CINT, cint)
    addmat!(:DIFFUSION, diffusion); addmat!(:LAMBDA, lambda); addmat!(:Jy, jy)
    addmat!(:MANIFESTMEANS, manifestmeans); addmat!(:MANIFESTVAR, manifestvar)
    addmat!(:T0VAR, t0var); addmat!(:T0MEANS, t0means); addmat!(:PARS, pars)

    DataFrame(matrix=matrices, row=rows, col=cols, parnumber=parnumber, value=value,
        transform=transform, predicttransform=predicttransform,
        updatetransform=fill(missing, length(matrices)))
end

"""Panel data: `nsubjects` subjects of slightly differing length."""
function _laplace_test_data(nsubjects, nmanifest; seed=11)
    subject_starts = Int[]; times = Float64[]; columns = Vector{Float64}[]
    position = 1
    for s in 1:nsubjects
        push!(subject_starts, position)
        nobs = 4 + (s % 3)
        for t in 1:nobs
            push!(times, 0.5 * (t - 1))
            push!(columns, [0.4 * sin(seed + s + 1.7t + 0.3m) for m in 1:nmanifest])
            position += 1
        end
    end
    return subject_starts, times, reduce(hcat, columns)
end

# A linear one-latent model whose T0MEANS and CINT are free, identity-transformed
# parameters -- so the log likelihood is exactly quadratic in any shift applied
# to them, which is what makes the Gaussian reference below available.
#
# Raw vector layout: 1-5 are the model parameters, 6-7 the population scales for
# the two random effects, 8 their correlation. The population parameters live in
# a tail of the same vector that no matrix cell reads, which is exactly how the
# R side lays them out.
function _build_laplace_linear_objective(; nsubjects=6)
    df = _laplace_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.4;;], t0means=[0.0;;],
        free=Dict(
            (:T0MEANS, 1, 1) => (1, "param[1]"),
            (:CINT, 1, 1) => (2, "param[2]"),
            (:DRIFT, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:JAx, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:MANIFESTMEANS, 1, 1) => (5, "param[5]"),
        ),
    )
    sp = ekf_from_data_frame(df)
    subject_starts, times, data = _laplace_test_data(nsubjects, 1)
    return ctsem_objective(sp, subject_starts, times, data)
end

# Each model *shape* is built exactly once and shared, for the reason
# test_adjoint_gradient_validation.jl gives: `ekf_from_data_frame` parses
# transform strings through `eval`, so every call produces uniquely-typed
# closures and forces the whole EKF pipeline to compile again for that type.
# Here that cost is paid three times over -- primal, dual, and the nested dual
# the exact outer gradient needs -- so rebuilding per testset dominated this
# file's runtime entirely. What each testset actually needs to be independent
# is a fresh *Laplace wrapper*, since the inner modes are the only state; that
# is free, and `_fresh_linear` provides it.
const _LAPLACE_LINEAR_OBJECTIVE = _build_laplace_linear_objective()
const _LAPLACE_LINEAR_VALUES = [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4]

"""A Laplace object with no retained modes, over the shared linear objective."""
_fresh_linear() = (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
    [1, 2], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))

# The same shape, but with a state-dependent DRIFT. The integrand is no longer
# Gaussian in `z`, so the Laplace value is genuinely approximate -- which is
# precisely why the gradient check has to pass here too.
function _build_laplace_nonlinear_objective(; nsubjects=5)
    expr = "PARS[1,1] * (1 + 0.15 * state[1])"
    df = _laplace_test_dataframe(
        drift=[0.0;;], jax=[0.0;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.4;;], t0means=[0.0;;], pars=[0.0;;],
        free=Dict(
            (:T0MEANS, 1, 1) => (1, "param[1]"),
            (:PARS, 1, 1) => (2, "-log1p_exp(param[2])"),
            (:DIFFUSION, 1, 1) => (3, "log1p_exp(param[3])"),
            (:MANIFESTMEANS, 1, 1) => (4, "param[4]"),
        ),
        predict=Dict((:DRIFT, 1, 1) => expr, (:JAx, 1, 1) => expr),
    )
    sp = ekf_from_data_frame(df)
    subject_starts, times, data = _laplace_test_data(nsubjects, 1; seed=5)
    return ctsem_objective(sp, subject_starts, times, data)
end

const _LAPLACE_NONLINEAR_OBJECTIVE = _build_laplace_nonlinear_objective()

# One random effect, on T0MEANS. Raw position 5 is its population scale.
_fresh_nonlinear() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE,
    [1], [5], Int[], [1.0]), [0.1, 0.4, -0.2, 0.05, -0.25])

################################################################################
# Independent references
################################################################################

"""
    _subject_loglik(laplace, i, values, z)

Subject `i`'s process log likelihood at the shifted parameter vector, through
the ordinary primal objective -- no adjoint, no tape, no Laplace machinery.
"""
function _subject_loglik(laplace, i, values, z)
    spec = laplace.spec
    L = ContinuousTimeSEM._laplace_popchol(collect(Float64, values), spec)
    shifted = ContinuousTimeSEM._laplace_subject_values(collect(Float64, values), spec, L, collect(Float64, z))
    return laplace.objective.subject_objectives[i](shifted)
end

"""
    _gaussian_reference(laplace, i, values)

`log integral exp(ll_i(z)) N(z|0,I) dz`, in closed form, for a model whose log
likelihood is exactly quadratic in `z`.

The quadratic is recovered from the log likelihood by evaluation -- value at
the origin, at each axis, and at each pair -- rather than assumed, so this
reference is derived from the same likelihood the Laplace code sees but by a
completely different route. `_assert_quadratic` checks the recovery is faithful
before the reference is trusted.

With `ll(z) = c + b'z + z'Az` and `P = I/2 - A`,

    log integral = c + b'P^-1 b / 4 - logdet(2P) / 2
"""
function _gaussian_reference(laplace, i, values)
    k = ContinuousTimeSEM.nrandomeffects(laplace.spec)
    ll(z) = _subject_loglik(laplace, i, values, z)
    c = ll(zeros(k))
    b = zeros(k)
    A = zeros(k, k)
    for j in 1:k
        e = zeros(k); e[j] = 1.0
        plus = ll(e); minus = ll(-e)
        b[j] = (plus - minus) / 2
        A[j, j] = (plus + minus) / 2 - c
    end
    for j in 1:k, l in (j + 1):k
        e = zeros(k); e[j] = 1.0; e[l] = 1.0
        off = ll(e) - c - b[j] - b[l] - A[j, j] - A[l, l]
        A[j, l] = off / 2
        A[l, j] = off / 2
    end
    P = Matrix(I / 2 - A)
    return (value=c + dot(b, P \ b) / 4 - logdet(2 .* P) / 2, c=c, b=b, A=A)
end

"""Confirm the recovered quadratic reproduces the log likelihood off its knots."""
function _assert_quadratic(laplace, i, values, reference)
    k = ContinuousTimeSEM.nrandomeffects(laplace.spec)
    for probe in (fill(0.7, k), fill(-1.3, k), collect(range(-0.4, 0.9, length=k)))
        predicted = reference.c + dot(reference.b, probe) + dot(probe, reference.A * probe)
        @test isapprox(predicted, _subject_loglik(laplace, i, values, probe); rtol=1e-9)
    end
end

"""Central difference of the Laplace *value*, the referee for its gradient."""
function _value_finite_difference(laplace, values; step=1e-5, gradient_method=:exact)
    x = collect(Float64, values)
    out = similar(x)
    for j in eachindex(x)
        h = step * max(1.0, abs(x[j]))
        plus = copy(x); plus[j] += h
        minus = copy(x); minus[j] -= h
        vp = ctsem_laplace_evaluate(laplace, plus; gradient=false, gradient_method=gradient_method).value
        vm = ctsem_laplace_evaluate(laplace, minus; gradient=false, gradient_method=gradient_method).value
        out[j] = (vp - vm) / (2h)
    end
    return out
end

################################################################################
# Tests
################################################################################

@testset "Laplace is exact for a Gaussian integrand" begin
    laplace, values = _fresh_linear()
    nsubjects = length(laplace.objective.subject_objectives)

    # The Laplace total, less the terms that are not per-subject integrals.
    result = ctsem_laplace_evaluate(laplace, values; gradient=false)
    @test isfinite(result.value)
    @test result.converged

    reference_total = 0.0
    for i in 1:nsubjects
        reference = _gaussian_reference(laplace, i, values)
        _assert_quadratic(laplace, i, values, reference)
        reference_total += reference.value
    end
    # No prior is attached to this objective, so the value is the sum of the
    # per-subject integrals and nothing else.
    @test isapprox(result.value, reference_total; rtol=1e-9)
end

@testset "the inner mode is a mode" begin
    laplace, values = _fresh_linear()
    ctsem_laplace_evaluate(laplace, values; gradient=false)
    diagnostics = ctsem_laplace_diagnostics(laplace)
    @test all(diagnostics.converged)
    @test all(diagnostics.max_gradient .< 1e-9)
    @test !any(diagnostics.hessian_repaired)

    # A directly checked stationarity condition: perturbing z away from the
    # reported mode may not increase the inner objective.
    spec = laplace.spec
    k = ContinuousTimeSEM.nrandomeffects(spec)
    for i in 1:length(laplace.objective.subject_objectives)
        zhat = Vector{Float64}(laplace.modes[:, i])
        at_mode = _subject_loglik(laplace, i, values, zhat) - dot(zhat, zhat) / 2
        for j in 1:k, delta in (-0.05, 0.05)
            probe = copy(zhat); probe[j] += delta
            nearby = _subject_loglik(laplace, i, values, probe) - dot(probe, probe) / 2
            @test nearby <= at_mode + 1e-10
        end
    end
end

@testset "the exact outer gradient is the gradient of the value" begin
    laplace, values = _fresh_linear()
    result = ctsem_laplace_evaluate(laplace, values; gradient=true, gradient_method=:exact)
    @test !result.approximate
    reference = _value_finite_difference(laplace, values)
    @test norm(result.gradient - reference) / norm(reference) < 1e-6
end

@testset "the exact gradient holds for a nonlinear model too" begin
    # The value here is genuinely approximate -- the integrand is not Gaussian
    # in z. The gradient still has to be the gradient *of that value*, which is
    # what makes the outer optimizer's convergence meaningful.
    laplace, values = _fresh_nonlinear()
    result = ctsem_laplace_evaluate(laplace, values; gradient=true, gradient_method=:exact)
    @test isfinite(result.value)
    reference = _value_finite_difference(laplace, values)
    @test norm(result.gradient - reference) / norm(reference) < 1e-5
end

@testset "the approximate gradient shares the value but drops the trace term" begin
    laplace, values = _fresh_linear()
    exact = ctsem_laplace_evaluate(laplace, values; gradient=true, gradient_method=:exact)
    approximate = ctsem_laplace_evaluate(laplace, values; gradient=true, gradient_method=:approximate)

    # Same objective, so the same value -- only the gradient is cheaper.
    @test isapprox(exact.value, approximate.value; rtol=1e-12)
    @test approximate.approximate
    @test !exact.approximate
    # And it is genuinely a different gradient, not a silently identical one:
    # if these agreed, the exact path would not be computing anything extra.
    @test norm(exact.gradient - approximate.gradient) / norm(exact.gradient) > 1e-6
end

@testset "the envelope gradient is the forward sweep it replaces" begin
    # The approximate gradient is assembled by hand from one reverse sweep per
    # subject plus the population-covariance Jacobian, rather than by
    # differentiating the summed inner objective in forward mode. The hand
    # assembly is what makes it cheap; this checks it is the same number.
    #
    # The referee holds the modes fixed, exactly as the envelope theorem
    # licenses, and differentiates the sum of inner objectives directly.
    laplace, values = _fresh_linear()
    result = ctsem_laplace_evaluate(laplace, values; gradient=true,
        gradient_method=:approximate)

    modes = [Vector{Float64}(laplace.modes[:, i])
             for i in 1:length(laplace.objective.subject_objectives)]
    summed_inner = function (x)
        S = eltype(x)
        ws = ContinuousTimeSEM._laplace_workspace!(laplace, S, length(x))
        Ld = ContinuousTimeSEM._laplace_popchol(x, laplace.spec)
        total = zero(S)
        for i in eachindex(modes)
            z = convert(Vector{S}, modes[i])
            total += ContinuousTimeSEM._laplace_inner_objective_gradient(
                laplace, i, x, Ld, z, ws).value
        end
        return total
    end
    reference = ForwardDiff.gradient(summed_inner, collect(values))
    @test norm(result.gradient - reference) / norm(reference) < 1e-9

    # And it is not accidentally the exact gradient: the dropped term is real.
    exact = ctsem_laplace_evaluate(laplace, values; gradient=true, gradient_method=:exact)
    @test norm(exact.gradient - result.gradient) / norm(exact.gradient) > 1e-6
end

@testset "the population covariance follows the Stan parameterisation" begin
    laplace, values = _fresh_linear()
    spec = laplace.spec
    L = ContinuousTimeSEM._laplace_popchol(collect(Float64, values), spec)
    @test size(L) == (2, 2)
    @test istril(L)

    # Reconstructed here the way ctModelWriter.R writes it, so a change to
    # either parameterisation shows up as a disagreement rather than as a
    # quietly different population distribution.
    scales = [ContinuousTimeSEM.log1p_exp(2 * values[spec.sd_index[j]] - 1) *
              spec.sd_scale[j] + 1e-10 for j in 1:2]
    base = zeros(2, 2)
    base[1, 1] = scales[1]; base[2, 2] = scales[2]
    base[2, 1] = 2 / (1 + exp(-values[spec.cor_index[1]])) - 1
    corsqrt = ContinuousTimeSEM.constraincorsqrt1(base)
    correlation = corsqrt * corsqrt'
    scaled = scales .+ 1e-8
    expected = (scaled .* correlation) .* scaled'
    @test isapprox(L * L', expected; rtol=1e-12)

    covariance = ctsem_laplace_popcov(laplace, values)
    @test isapprox(covariance, L * L'; rtol=1e-12)
    @test isposdef(Symmetric(covariance))
end

@testset "per-subject modes and conditional standard errors" begin
    laplace, values = _fresh_linear()
    nsubjects = length(laplace.objective.subject_objectives)
    modes = ctsem_laplace_modes(laplace, values)

    @test size(modes.z) == (nsubjects, 2)
    @test size(modes.raw) == (nsubjects, 2)
    @test all(modes.converged)
    @test all(modes.z_sd .> 0)
    @test all(modes.raw_sd .> 0)

    # The raw-scale mode is the standardized one pushed through L, which is the
    # same change of variables the model applies to produce a subject's
    # parameters.
    L = ContinuousTimeSEM._laplace_popchol(collect(Float64, values), laplace.spec)
    for i in 1:nsubjects
        @test isapprox(modes.raw[i, :], L * modes.z[i, :]; rtol=1e-10)
    end
    # Subjects differ: a collapsed mode vector would mean the random effects
    # are not being identified at all.
    @test maximum(abs, diff(modes.z; dims=1)) > 1e-6
end

@testset "warm starting does not change the answer" begin
    # The modes persist between calls, so a second evaluation at the same point
    # starts from the answer rather than from zero. That must be an efficiency,
    # not a difference.
    laplace, values = _fresh_linear()
    first_result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    second_result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    @test isapprox(first_result.value, second_result.value; rtol=1e-12)
    @test isapprox(first_result.gradient, second_result.gradient; rtol=1e-10)

    fresh, _ = _fresh_linear()
    fresh_result = ctsem_laplace_evaluate(fresh, values; gradient=true)
    @test isapprox(fresh_result.value, first_result.value; rtol=1e-10)
end

@testset "zero random effects reduces to the ordinary likelihood" begin
    # The degenerate case has to be a no-op rather than an error: it is what a
    # model with the Laplace path selected but nothing declared varying is.
    df = _laplace_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.4;;], t0means=[0.0;;],
        free=Dict(
            (:T0MEANS, 1, 1) => (1, "param[1]"),
            (:DRIFT, 1, 1) => (2, "-log1p_exp(param[2])"),
            (:JAx, 1, 1) => (2, "-log1p_exp(param[2])"),
        ),
    )
    sp = ekf_from_data_frame(df)
    subject_starts, times, data = _laplace_test_data(4, 1)
    objective = ctsem_objective(sp, subject_starts, times, data)
    laplace = ctsem_laplace_objective(objective, Int[], Int[], Int[], Float64[])
    values = [0.2, 0.3]

    result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    @test isapprox(result.value, objective(values); rtol=1e-12)
    @test isapprox(result.gradient, ctsem_adjoint_gradient(objective, values).gradient; rtol=1e-8)
end

@testset "the outer Hessian is symmetric and matches the gradient" begin
    laplace, values = _fresh_linear()
    H = ctsem_laplace_hessian(laplace, values; step=1e-4)
    n = length(values)
    @test size(H) == (n, n)
    @test H == transpose(H)
    @test all(isfinite, H)

    # Referee: a directional derivative of the exact gradient.
    direction = normalize(collect(1.0:n))
    h = 1e-5
    gp = ctsem_laplace_evaluate(laplace, values .+ h .* direction; gradient=true).gradient
    gm = ctsem_laplace_evaluate(laplace, values .- h .* direction; gradient=true).gradient
    @test norm(H * direction - (gp .- gm) ./ (2h)) / norm(H * direction) < 1e-4
end

@testset "an unusable inner curvature is reported, not hidden" begin
    laplace, values = _fresh_linear()
    # A population scale driven to essentially zero makes the random effects
    # unidentified; the run must still produce a finite value and say what
    # happened, rather than returning a silently meaningless number.
    degenerate = copy(values)
    degenerate[laplace.spec.sd_index] .= -40.0
    result = ctsem_laplace_evaluate(laplace, degenerate; gradient=false)
    @test isfinite(result.value)
    diagnostics = ctsem_laplace_diagnostics(laplace)
    @test length(diagnostics.converged) == length(laplace.objective.subject_objectives)
    @test all(isfinite, diagnostics.max_gradient)
end
