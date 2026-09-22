# Model fixtures shared by the adjoint tests and `test_discrete_time.jl`.
#
# In a file of their own because the parallel runner gives each test file its
# own process while `runtests.jl` gives them one: a fixture defined in whichever
# file happens to be included first is there under one runner and missing under
# the other, and it reports as an `UndefVarError` naming the fixture rather than
# as a missing include.

function _adjoint_test_dataframe(;
    drift::AbstractMatrix, jax::AbstractMatrix, cint::AbstractMatrix,
    diffusion::AbstractMatrix, lambda::AbstractMatrix, jy::AbstractMatrix,
    manifestmeans::AbstractMatrix, manifestvar::AbstractMatrix,
    t0var::AbstractMatrix, t0means::AbstractMatrix,
    pars::AbstractMatrix=zeros(1, 1),
    tdpredeffect::Union{Nothing,AbstractMatrix}=nothing,
    jtd::Union{Nothing,AbstractMatrix}=nothing,
    free::Dict=Dict(), predict::Dict=Dict(), update::Dict=Dict(),
)
    matrices = Symbol[]
    rows = Int[]
    cols = Int[]
    parnumber = Union{Missing,Int}[]
    value = Union{Missing,Float64}[]
    transform = Union{Missing,String}[]
    predicttransform = Union{Missing,String}[]
    updatetransform = Union{Missing,String}[]

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
            push!(updatetransform, get(update, key, missing))
        end
    end

    addmat!(:DRIFT, drift); addmat!(:JAx, jax); addmat!(:CINT, cint)
    addmat!(:DIFFUSION, diffusion); addmat!(:LAMBDA, lambda); addmat!(:Jy, jy)
    addmat!(:MANIFESTMEANS, manifestmeans); addmat!(:MANIFESTVAR, manifestvar)
    addmat!(:T0VAR, t0var); addmat!(:T0MEANS, t0means); addmat!(:PARS, pars)
    if tdpredeffect !== nothing
        addmat!(:TDPREDEFFECT, tdpredeffect); addmat!(:Jtd, jtd)
    end

    DataFrame(matrix=matrices, row=rows, col=cols, parnumber=parnumber, value=value,
        transform=transform, predicttransform=predicttransform, updatetransform=updatetransform)
end

# Ties DRIFT[1,1] and JAx[1,1] to one negative-constrained free parameter and
# DIFFUSION[1,1] to one positive-constrained free parameter -- the smallest
# genuinely free (not all-fixed) linear model.
function _adjoint_linear_1d_parameters()
    df = _adjoint_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DIFFUSION, 1, 1) => (2, "log1p_exp(param[2])"),
        ),
    )
    ekf_from_data_frame(df)
end

# A 2-latent, 2-manifest model with a free cross-effect, so the gradient
# check exercises coupled dynamics, not just a diagonal system.
#
# `continuous` is threaded through to `ekf_from_data_frame` so
# `test_discrete_time.jl` can reuse this exact cell layout for its discrete
# variant (F5) instead of writing a new fixture.
function _adjoint_cross_effect_2d_parameters(; continuous::Bool=true)
    df = _adjoint_test_dataframe(
        drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3],
        cint=[0.0; 0.0;;], diffusion=[0.2 0.0; 0.0 0.15],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.0 0.1],
        t0var=[1.0 0.0; 0.0 1.0], t0means=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DRIFT, 2, 2) => (2, "-log1p_exp(param[2])"),
            (:JAx, 2, 2) => (2, "-log1p_exp(param[2])"),
            (:DRIFT, 1, 2) => (3, "param[3]"),
            (:JAx, 1, 2) => (3, "param[3]"),
        ),
    )
    ekf_from_data_frame(df; continuous_time=continuous)
end

# PARS[1,1] is a free parameter feeding a state-dependent DRIFT/JAx
# expression -- exercises `apply_complex_transforms_at_indices!` under
# ForwardDiff, which is exactly the code path a hand-written or Enzyme-based
# adjoint has to reproduce correctly for nonlinear ctsem models.
function _adjoint_state_dependent_1d_parameters()
    expr = "PARS[1,1] * (1 + 0.1 * state[1])"
    df = _adjoint_test_dataframe(
        drift=[0.0;;], jax=[0.0;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.5;;], t0means=[0.0;;], pars=[0.0;;],
        free=Dict((:PARS, 1, 1) => (1, "param[1]")),
        predict=Dict((:DRIFT, 1, 1) => expr, (:JAx, 1, 1) => expr),
    )
    ekf_from_data_frame(df)
end

# Each of these three model *types* is built exactly once and reused across
# every @testset that needs it. `ekf_from_data_frame` parses transform
# strings via `eval`, so every call produces uniquely-typed closures and
# forces a full fresh compilation of the EKF pipeline for that type -- fine
# for real model fitting (once per R session) but wasteful if paid once per
# @testset here for what is really the same 3 model types varying only in
# data. Building each sp once cuts this file's compilation cost roughly in
# half versus one fresh `_adjoint_*_parameters()` call per testset.
# Everything the six original scenarios leave fixed -- T0VAR, MANIFESTVAR,
# MANIFESTMEANS, CINT, and the off-diagonal correlation-sqrt entries -- is free
# here, so the reverse paths through `sdcovsqrt2cov` (for all three covariance
# matrices), the discrete-intercept solve, and the measurement-mean term are
# actually exercised rather than multiplied by a structurally zero cotangent.
#
# `continuous` is threaded through for the same reason as in
# `_adjoint_cross_effect_2d_parameters` above (F5).
#
# F6: LAMBDA[2,1] is also free (parameter 11), so `_reverse_update!`'s
# `Λ̄` scatter (`adjoint_ekf.jl:842`) is multiplied by something. Before this,
# LAMBDA and Jy carried `parnumber = missing` in every Julia test model, so a
# transposition, sign or orientation error in the LAMBDA/Jy/Jtd cotangent
# scatters was multiplied by nothing across the whole suite.
function _adjoint_free_covariance_2d_parameters(; continuous::Bool=true)
    df = _adjoint_test_dataframe(
        drift=[-0.5 0.3; 0.1 -0.3], jax=[-0.5 0.3; 0.1 -0.3],
        cint=[0.0; 0.0;;], diffusion=[0.2 0.0; 0.05 0.15],
        lambda=[1.0 0.0; 0.0 1.0], jy=[1.0 0.0; 0.0 1.0],
        manifestmeans=[0.0; 0.0;;], manifestvar=[0.1 0.0; 0.02 0.1],
        t0var=[1.0 0.0; 0.03 1.0], t0means=[0.0; 0.0;;],
        free=Dict(
            (:DRIFT, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:JAx, 1, 1) => (1, "-log1p_exp(param[1])"),
            (:DRIFT, 1, 2) => (2, "param[2]"),
            (:JAx, 1, 2) => (2, "param[2]"),
            (:CINT, 1, 1) => (3, "param[3]"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:DIFFUSION, 2, 1) => (5, "param[5]"),
            (:MANIFESTVAR, 2, 2) => (6, "log1p_exp(param[6])"),
            (:MANIFESTMEANS, 1, 1) => (7, "param[7]"),
            (:T0VAR, 1, 1) => (8, "log1p_exp(param[8])"),
            (:T0VAR, 2, 1) => (9, "param[9]"),
            (:T0MEANS, 2, 1) => (10, "param[10]"),
            (:LAMBDA, 2, 1) => (11, "param[11]"),
        ),
    )
    ekf_from_data_frame(df; continuous_time=continuous)
end
