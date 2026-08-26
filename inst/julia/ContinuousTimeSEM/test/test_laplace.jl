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
function _value_finite_difference(laplace, values; step=1e-5, )
    x = collect(Float64, values)
    out = similar(x)
    for j in eachindex(x)
        h = step * max(1.0, abs(x[j]))
        plus = copy(x); plus[j] += h
        minus = copy(x); minus[j] -= h
        vp = ctsem_laplace_evaluate(laplace, plus; gradient=false).value
        vm = ctsem_laplace_evaluate(laplace, minus; gradient=false).value
        out[j] = (vp - vm) / (2h)
    end
    return out
end


# --- two levels: subjects nested in studies ----------------------------------
#
# Six subjects in three studies. The subject level carries T0MEANS and CINT,
# the study level carries MANIFESTMEANS -- all identity transforms on the
# shared linear objective, so the integrand is exactly Gaussian in the unit
# latent vector and Laplace is exact here too, which is what makes the
# closed-form reference below available for the *coupled* case.
#
# Raw layout: 1-5 the model parameters, 6-7 the subject scales, 8 their
# correlation, 9 the study scale.
_TWOLEVEL_GROUP = [1, 1, 2, 2, 3, 3]
_fresh_twolevel() = (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE;
    re_index=[1, 2, 5], sd_index=[6, 7, 9], cor_index=[8],
    sd_scale=[1.0, 1.0, 1.0], level_nre=[2, 1],
    group=vcat(1:6, _TWOLEVEL_GROUP), level_ngroups=[6, 3]),
    [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25])


# --- three levels: subjects in studies in regions ----------------------------
#
# The same six subjects, their three studies grouped into two regions. Nothing
# about the assembly is specific to two levels -- a block simply acquires one
# cross term per ancestor -- but "one ancestor" and "a chain of ancestors" are
# different code paths through the selected inverse and through the explicit
# `dL` terms, so the depth that exercises both is tested rather than assumed.
#
# Raw layout: 1-5 the model parameters, 6-7 the subject scales, 8 their
# correlation, 9 the study scale, 10 the region scale.
_THREELEVEL_STUDY = [1, 1, 2, 2, 3, 3]
_THREELEVEL_REGION = [1, 1, 1, 1, 2, 2]
_fresh_threelevel() = (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE;
    re_index=[1, 2, 5, 3], sd_index=[6, 7, 9, 10], cor_index=[8],
    sd_scale=[1.0, 1.0, 1.0, 1.0], level_nre=[2, 1, 1],
    group=vcat(1:6, _THREELEVEL_STUDY, _THREELEVEL_REGION),
    level_ngroups=[6, 3, 2]),
    [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25, 0.1])

"""Unit `U`'s log likelihood at latent vector `u`, through the primal only."""
function _unit_loglik(laplace, U, values, u)
    spec = laplace.spec
    Ls = ContinuousTimeSEM._laplace_popchols(collect(Float64, values), spec)
    total = 0.0
    for (m, i) in enumerate(laplace.units.members[U])
        shifted = ContinuousTimeSEM._laplace_member_values(collect(Float64, values),
            spec, Ls, collect(Float64, u), laplace.units.offsets[U][m])
        total += laplace.objective.subject_objectives[i](shifted)
    end
    return total
end

"""
Closed-form `log integral exp(ll_U(u)) N(u|0,I) du` for a unit whose log
likelihood is exactly quadratic in `u`, recovered from the log likelihood by
evaluation rather than assumed.
"""
function _unit_gaussian_reference(laplace, U, values)
    d = laplace.units.dims[U]
    ll(u) = _unit_loglik(laplace, U, values, u)
    c = ll(zeros(d))
    b = zeros(d); A = zeros(d, d)
    for j in 1:d
        e = zeros(d); e[j] = 1.0
        plus = ll(e); minus = ll(-e)
        b[j] = (plus - minus) / 2
        A[j, j] = (plus + minus) / 2 - c
    end
    for j in 1:d, l in (j + 1):d
        e = zeros(d); e[j] = 1.0; e[l] = 1.0
        off = ll(e) - c - b[j] - b[l] - A[j, j] - A[l, l]
        A[j, l] = off / 2; A[l, j] = off / 2
    end
    P = Matrix(I / 2 - A)
    return c + dot(b, P \ b) / 4 - logdet(2 .* P) / 2
end

################################################################################
# Tests
################################################################################

@testset "two levels: units couple the subjects that share a study" begin
    laplace, values = _fresh_twolevel()
    units = laplace.units
    # Three studies of two subjects each, so three units rather than six.
    @test length(units.members) == 3
    @test sort(vcat(units.members...)) == collect(1:6)
    for U in 1:3
        @test length(units.members[U]) == 2
        # Two subject blocks of 2 plus one shared study block of 1.
        @test units.dims[U] == 2 * 2 + 1
        # Both members point at the *same* study block -- that sharing is the
        # coupling, and it is the thing that makes a study one integration unit.
        @test units.offsets[U][1][2] == units.offsets[U][2][2]
        # ...and at different subject blocks.
        @test units.offsets[U][1][1] != units.offsets[U][2][1]
    end
end

@testset "two levels: Laplace is exact for a Gaussian integrand" begin
    # The same argument as the single-level exactness test, but over a unit
    # latent vector that mixes two subject blocks and a shared study block. If
    # the offsets, the per-level Cholesky factors or the shared block were
    # wrong, this is where it shows.
    laplace, values = _fresh_twolevel()
    result = ctsem_laplace_evaluate(laplace, values; gradient=false)
    @test isfinite(result.value)
    @test result.converged

    reference = sum(_unit_gaussian_reference(laplace, U, values) for U in 1:3)
    @test isapprox(result.value, reference; rtol=1e-8)
end

@testset "two levels: the outer gradient is the gradient of the value" begin
    laplace, values = _fresh_twolevel()
    result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    reference = _value_finite_difference(laplace, values)
    @test norm(result.gradient - reference) / norm(reference) < 1e-6
end

@testset "two levels: an empty outer level reduces to one level" begin
    # A study level carrying no random effects is a hierarchy that does not do
    # anything, and it must not change the answer -- only how the subjects are
    # grouped for integration. This is the degenerate case that catches an
    # offset or dimension bug in the unit layout.
    single, values = _fresh_linear()
    grouped = ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE;
        re_index=[1, 2], sd_index=[6, 7], cor_index=[8], sd_scale=[1.0, 1.0],
        level_nre=[2, 0], group=vcat(1:6, _TWOLEVEL_GROUP), level_ngroups=[6, 3])

    a = ctsem_laplace_evaluate(single, values; gradient=true, nested_gradient=true)
    b = ctsem_laplace_evaluate(grouped, values; gradient=true)
    @test isapprox(a.value, b.value; rtol=1e-9)
    @test norm(a.gradient - b.gradient) / norm(a.gradient) < 1e-7
end

@testset "two levels: modes are reported per level and per group" begin
    laplace, values = _fresh_twolevel()
    ctsem_laplace_evaluate(laplace, values; gradient=false)

    subject = ctsem_laplace_modes(laplace, values, 1)
    study = ctsem_laplace_modes(laplace, values, 2)
    # One row per subject at level 1, one per *study* at level 2 -- a study
    # effect is one vector shared by its members, not one per member.
    @test size(subject.z) == (6, 2)
    @test size(study.z) == (3, 1)
    @test study.ngroups == 3
    @test study.group == _TWOLEVEL_GROUP
    @test all(isfinite, subject.z)
    @test all(isfinite, study.z)
    @test all(study.z_sd .> 0)
    # Studies differ; a collapsed study mode would mean the level is doing
    # nothing.
    @test maximum(abs, diff(study.z; dims=1)) > 1e-8
end

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
        zhat = laplace.modes[i]
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
    result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    reference = _value_finite_difference(laplace, values)
    @test norm(result.gradient - reference) / norm(reference) < 1e-6
end

@testset "the exact gradient holds for a nonlinear model too" begin
    # The value here is genuinely approximate -- the integrand is not Gaussian
    # in z. The gradient still has to be the gradient *of that value*, which is
    # what makes the outer optimizer's convergence meaningful.
    laplace, values = _fresh_nonlinear()
    result = ctsem_laplace_evaluate(laplace, values; gradient=true)
    @test isfinite(result.value)
    reference = _value_finite_difference(laplace, values)
    @test norm(result.gradient - reference) / norm(reference) < 1e-5
end

@testset "the seeded gradient matches the nested one it replaces" begin
    # The production gradient assembles `dT/dtheta` from `k + 1` seeded reverse
    # sweeps plus a good deal of chain rule. The nested route computes the same
    # quantity by running ForwardDiff over the entire per-subject term. They
    # share the primal and nothing else -- no seeding, no hand-written
    # assembly, no popchol derivatives -- so agreement to machine precision is
    # a real check on every term of that assembly, and much sharper than the
    # finite-difference comparison can be.
    for (label, fresh) in (("linear", _fresh_linear), ("nonlinear", _fresh_nonlinear))
        laplace, values = fresh()
        seeded = ctsem_laplace_evaluate(laplace, values; gradient=true)
        nested = ctsem_laplace_evaluate(laplace, values; gradient=true,
            nested_gradient=true)
        @test seeded.value == nested.value
        @test norm(seeded.gradient - nested.gradient) / norm(nested.gradient) < 1e-8
    end
end

@testset "the seeded gradient is exercised away from the optimum" begin
    # A single test point can hide a term that happens to vanish there. These
    # perturb the population scales and correlations, which are the parameters
    # the assembly treats specially -- they move `v` through `L` as well as
    # directly, and they are the only ones with an explicit `psi` term.
    laplace, values = _fresh_linear()
    spec = laplace.spec
    for shift in (0.4, -0.6, 1.1)
        probe = collect(values)
        probe[spec.levels[1].sd_index] .+= shift
        probe[spec.levels[1].cor_index] .-= shift / 2
        seeded = ctsem_laplace_evaluate(laplace, probe; gradient=true)
        nested = ctsem_laplace_evaluate(laplace, probe; gradient=true,
            nested_gradient=true)
        @test norm(seeded.gradient - nested.gradient) / norm(nested.gradient) < 1e-8
    end
end

@testset "chunking the subject loop does not change the answer" begin
    # `cores` reaches this path as a chunk count, and each chunk owns its own
    # adjoint workspace rather than indexing one by thread id -- a task can
    # migrate between threads at any yield point, so thread-indexed scratch
    # would be a race. Splitting only reorders a sum, so the tolerance is
    # floating-point association, not an approximation.
    laplace, values = _fresh_linear()
    original = ctsem_max_chunks().max_chunks
    try
        ctsem_set_max_chunks!(1)
        serial = ctsem_laplace_evaluate(laplace, values; gradient=true)
        for chunks in (2, 3, 8)
            fresh, _ = _fresh_linear()
            ctsem_set_max_chunks!(chunks)
            split = ctsem_laplace_evaluate(fresh, values; gradient=true)
            @test isapprox(split.value, serial.value; rtol=1e-12)
            @test isapprox(split.gradient, serial.gradient; rtol=1e-9)
            @test isapprox(split.subject_loglik, serial.subject_loglik; rtol=1e-12)
        end
    finally
        ctsem_set_max_chunks!(original)
    end
end

@testset "the blocked curvature equals the dense one" begin
    # The production curvature is assembled one block of `u` at a time, using
    # the fact that a member outside a block has structurally zero second
    # derivative with respect to it. The dense route differentiates the whole
    # gradient at once and knows nothing about that structure, so agreement is
    # a real check on the block bookkeeping -- which members own which block,
    # and where each block sits in `u`.
    for (label, fresh) in (("one level", _fresh_linear), ("two levels", _fresh_twolevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        Ls = ContinuousTimeSEM._laplace_popchols(theta, laplace.spec)
        for U in eachindex(laplace.units.members)
            u = 0.1 .* collect(1.0:laplace.units.dims[U])
            blocked = ContinuousTimeSEM._laplace_unit_hessian(laplace, U, theta, Ls, u)
            reference = ContinuousTimeSEM._laplace_unit_hessian(laplace, U, theta, Ls, u;
                dense=true)
            @test size(blocked) == size(reference)
            @test norm(blocked - reference) / max(norm(reference), 1) < 1e-10
        end
    end
end

@testset "the curvature is arrow structured, and the blocks say where" begin
    # Subjects in a study are conditionally independent given the study effect,
    # so their blocks do not couple to each other -- only to the study block.
    # That is the sparsity the blocked assembly exploits, so it is worth
    # asserting rather than assuming.
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    Ls = ContinuousTimeSEM._laplace_popchols(theta, laplace.spec)
    U = 1
    u = 0.1 .* collect(1.0:laplace.units.dims[U])
    H = ContinuousTimeSEM._laplace_unit_hessian(laplace, U, theta, Ls, u)

    blocks = laplace.units.blocks[U]
    subjectblocks = [b for b in blocks if length(b.members) == 1]
    studyblocks = [b for b in blocks if length(b.members) > 1]
    @test length(subjectblocks) == 2      # two subjects per study
    @test length(studyblocks) == 1        # one shared study block

    for a in subjectblocks, b in subjectblocks
        a === b && continue
        rows = (a.offset + 1):(a.offset + a.size)
        cols = (b.offset + 1):(b.offset + b.size)
        @test all(abs.(H[rows, cols]) .< 1e-12)
    end
    # ...and the coupling to the study block is genuinely there, so the test
    # above is not passing because everything is zero.
    sb = studyblocks[1]
    srows = (sb.offset + 1):(sb.offset + sb.size)
    a = subjectblocks[1]
    arows = (a.offset + 1):(a.offset + a.size)
    @test maximum(abs, H[arows, srows]) > 1e-8
end

@testset "correlations are capped at the boundary, and the cap is reported" begin
    # Driving the coordinate far past the cap must give the same covariance as
    # sitting exactly on it: that is what "capped" means, and it is what stops
    # the covariance running to singularity on a design that cannot identify
    # the correlation.
    laplace, values = _fresh_linear()
    cap = ContinuousTimeSEM._LAPLACE_COR_CAP[]
    spec = laplace.spec.levels[1]

    at = collect(values); at[spec.cor_index[1]] = cap
    past = collect(values); past[spec.cor_index[1]] = cap * 20
    @test ContinuousTimeSEM._laplace_popchol(at, spec) ≈
          ContinuousTimeSEM._laplace_popchol(past, spec)

    # And it is reported rather than silently applied.
    @test isempty(ctsem_laplace_boundary(laplace, values).level)
    flagged = ctsem_laplace_boundary(laplace, past)
    @test flagged.level == [1]
    @test flagged.position == [1]

    # Below the cap nothing is touched.
    below = collect(values); below[spec.cor_index[1]] = cap / 2
    @test !(ContinuousTimeSEM._laplace_popchol(below, spec) ≈
            ContinuousTimeSEM._laplace_popchol(at, spec))
    @test isempty(ctsem_laplace_boundary(laplace, below).level)
end

@testset "the mode Jacobian is the derivative it claims to be" begin
    # Subject parameters at posterior draws linearise the mode rather than
    # re-solving it. That is only defensible if the Jacobian is right and the
    # error is second order, so both are checked rather than asserted: halving
    # the step should quarter the discrepancy against a re-solved mode.
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    J = ctsem_laplace_mode_jacobian(laplace, theta)
    base = [copy(laplace.modes[U]) for U in eachindex(laplace.units.members)]

    direction = normalize(collect(1.0:length(theta)))
    errors = Float64[]
    for step in (0.02, 0.01)
        probe = theta .+ step .* direction
        Ls = ContinuousTimeSEM._laplace_popchols(probe, laplace.spec)
        worst = 0.0
        for U in eachindex(laplace.units.members)
            linear = base[U] .+ J[U] * (probe .- theta)
            # Re-solve from scratch at the probe, not warm-started from the
            # linearised guess, so the two are genuinely independent.
            laplace.modes[U] = zeros(length(base[U]))
            ContinuousTimeSEM._laplace_solve_unit_mode!(laplace, U, probe, Ls)
            worst = max(worst, maximum(abs, linear .- laplace.modes[U]))
        end
        push!(errors, worst)
        for U in eachindex(base); laplace.modes[U] = copy(base[U]); end
    end
    @test errors[1] < 1e-3
    # Second order: halving the step should cut the error by roughly four.
    @test errors[2] < errors[1] / 3
end

@testset "the block factorization equals the dense one" begin
    # The elimination exploits that a block couples only to its ancestors, so
    # it must reproduce what a dense Cholesky of the same matrix gives -- both
    # the log determinant and the solve. The dense route knows nothing about
    # the block tree, which is what makes agreement meaningful.
    laplace, values = _fresh_twolevel()
    theta = collect(Float64, values)
    Ls = ContinuousTimeSEM._laplace_popchols(theta, laplace.spec)
    for U in eachindex(laplace.units.members)
        blocks = laplace.units.blocks[U]
        u = 0.1 .* collect(1.0:laplace.units.dims[U])
        H = ContinuousTimeSEM._laplace_unit_hessian(laplace, U, theta, Ls, u)
        _, M = ContinuousTimeSEM._laplace_negate_definite(H)

        blocked = ContinuousTimeSEM._laplace_block_of(M, blocks)
        # Round tripping through the block form must not lose anything: if it
        # did, the sparsity pattern would be wrong rather than the arithmetic.
        @test ContinuousTimeSEM._laplace_block_dense(blocked, blocks,
            size(M, 1)) ≈ M

        ok, ld, factors, coupling = ContinuousTimeSEM._laplace_block_factor(blocked, blocks)
        @test ok
        @test ld ≈ logdet(cholesky(Symmetric(M)))

        rhs = collect(1.0:size(M, 1)) ./ size(M, 1)
        x = ContinuousTimeSEM._laplace_block_solve(factors, coupling, blocks, rhs)
        @test norm(M * x - rhs) / norm(rhs) < 1e-9
    end
end

@testset "the block factorization reports a curvature it cannot factor" begin
    laplace, values = _fresh_twolevel()
    blocks = laplace.units.blocks[1]
    M = ContinuousTimeSEM.CTSEMBlockMatrix(blocks)
    for d in M.diag; d .= -Matrix(I, size(d)...); end   # negative definite
    ok, ld, _, _ = ContinuousTimeSEM._laplace_block_factor(M, blocks)
    @test !ok
    @test isnan(ld)
end

@testset "the selected inverse matches the dense inverse on the pattern" begin
    # Only the entries inside `M`'s own sparsity pattern are produced, because
    # the full inverse of an arrow matrix is dense and is exactly what must not
    # be formed. Those entries have to be right, so they are compared with the
    # dense inverse of the same matrix.
    for (label, fresh) in (("one level", _fresh_linear), ("two levels", _fresh_twolevel))
        laplace, values = fresh()
        theta = collect(Float64, values)
        Ls = ContinuousTimeSEM._laplace_popchols(theta, laplace.spec)
        for U in eachindex(laplace.units.members)
            blocks = laplace.units.blocks[U]
            isempty(blocks) && continue
            u = 0.1 .* collect(1.0:laplace.units.dims[U])
            M = ContinuousTimeSEM._laplace_unit_curvature(laplace, U, theta, Ls, u)
            ContinuousTimeSEM._laplace_repair_blocks!(M, blocks)
            ok, _, factors, elim = ContinuousTimeSEM._laplace_block_factor(M, blocks)
            @test ok

            dense = ContinuousTimeSEM._laplace_block_dense(M, blocks,
                laplace.units.dims[U])
            reference = inv(Symmetric(dense))
            Cd, Cc = ContinuousTimeSEM._laplace_selected_inverse(factors, elim, blocks)

            for (b, block) in enumerate(blocks)
                rows = (block.offset + 1):(block.offset + block.size)
                @test norm(Cd[b] - reference[rows, rows]) /
                      max(norm(reference[rows, rows]), 1) < 1e-8
                for (t, a) in enumerate(block.ancestors)
                    cols = (blocks[a].offset + 1):(blocks[a].offset + blocks[a].size)
                    @test norm(Cc[b][t] - reference[rows, cols]) /
                          max(norm(reference[rows, cols]), 1) < 1e-8
                end
            end
        end
    end
end

@testset "the seeded nested gradient matches the nested oracle" begin
    # The production route for a hierarchy assembles `dT/dtheta` from
    # `O(members)` seeded sweeps plus the selected inverse. The nested route
    # computes the same quantity by running ForwardDiff over the whole
    # per-unit term, at a cost that scales with the parameter count. They
    # share the primal and nothing else, so agreement pins every term of the
    # assembly.
    for (label, fresh) in (("two levels", _fresh_twolevel),
                           ("three levels", _fresh_threelevel))
    laplace, values = fresh()
    theta = collect(Float64, values)
    Ls = ContinuousTimeSEM._laplace_popchols(theta, laplace.spec)
    dL = ContinuousTimeSEM._laplace_level_chol_derivatives(theta, laplace.spec)

    nunits = length(laplace.units.members)
    Ms = Vector{Any}(undef, nunits)
    curv = Vector{Any}(undef, nunits)
    for U in 1:nunits
        ContinuousTimeSEM._laplace_solve_unit_mode!(laplace, U, theta, Ls)
        M = ContinuousTimeSEM._laplace_unit_curvature(laplace, U, theta, Ls,
            laplace.modes[U])
        ContinuousTimeSEM._laplace_repair_blocks!(M, laplace.units.blocks[U])
        ok, _, f, e = ContinuousTimeSEM._laplace_block_factor(M, laplace.units.blocks[U])
        @test ok
        Ms[U] = M; curv[U] = (f, e)
    end

    seeded = zeros(Float64, length(theta))
    allok = true
    for U in 1:nunits
        allok &= ContinuousTimeSEM._laplace_seeded_unit_gradient!(seeded, laplace, U,
            theta, Ls, dL, Ms[U], curv[U][1], curv[U][2])
    end
    @test allok
    ContinuousTimeSEM._ctsem_log_prior_gradient!(seeded, laplace.objective, theta)

    nested = ContinuousTimeSEM._laplace_nested_gradient(laplace, theta, Ls, curv)
    @test (label, norm(seeded - nested) / norm(nested) < 1e-8) == (label, true)
    end
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
    scales = [ContinuousTimeSEM.log1p_exp(2 * values[spec.levels[1].sd_index[j]] - 1) *
              spec.levels[1].sd_scale[j] + 1e-10 for j in 1:2]
    base = zeros(2, 2)
    base[1, 1] = scales[1]; base[2, 2] = scales[2]
    base[2, 1] = 2 / (1 + exp(-values[spec.levels[1].cor_index[1]])) - 1
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
    degenerate[laplace.spec.levels[1].sd_index] .= -40.0
    result = ctsem_laplace_evaluate(laplace, degenerate; gradient=false)
    @test isfinite(result.value)
    diagnostics = ctsem_laplace_diagnostics(laplace)
    @test length(diagnostics.converged) == length(laplace.objective.subject_objectives)
    @test all(isfinite, diagnostics.max_gradient)
end
