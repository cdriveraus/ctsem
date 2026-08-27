# Model fixtures shared by the Laplace and quadrature suites.
#
# These were defined inside `test_laplace.jl` and reached from
# `test_quadrature.jl` by relying on `runtests.jl` including the files in sorted
# order into one module. That worked and was invisible until the suite started
# running one process per file, at which point the quadrature tests could no
# longer see them.
#
# So they live here, and both suites include this file. The `isdefined` guard is
# what lets the serial runner include it twice without redefining a `const`.
#
# Each model *shape* is built exactly once and shared, for the reason
# `test_adjoint_gradient_validation.jl` gives: `ekf_from_columns` parses
# transform strings through `eval`, so every call produces uniquely-typed
# closures and forces the whole EKF pipeline to compile again for that type.
# What a testset actually needs to be independent is a fresh *Laplace wrapper*,
# since the inner modes are the only state; that is free, and `_fresh_*`
# provides it.

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
