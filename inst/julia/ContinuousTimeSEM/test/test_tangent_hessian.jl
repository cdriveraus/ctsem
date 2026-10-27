using DataFrames, ForwardDiff, LinearAlgebra

# The hand-coded second order (`tangent_hessian.jl`) against the Hessian's dual
# route, its oracle: the two compute the same derivative of the same filter, so
# they agree to rounding. One model exercises what the hand route covers --
# every drift cell, a cross-loading, correlated diffusion and measurement
# errors, missing values and a row with none observed, substeps, unequal and
# zero subject weights, a prior, and two chunks of columns.

isdefined(@__MODULE__, :_adjoint_test_dataframe) ||
    include(joinpath(@__DIR__, "adjoint_fixtures.jl"))

function _tangent_test_model(; driftscale=1.0, max_timestep=Inf, perrow=false)
    drift = [-0.6 0.2; 0.1 -0.5]
    lambda = [1.0 0.0; 0.7 0.0; 0.0 1.0]
    # `driftscale` shrinks every drift cell: near zero, the filter integrates
    # the intervals as series.
    d(e) = driftscale == 1 ? e : "$(driftscale) * ($(e))"
    free = Dict(
        (:DRIFT, 1, 1) => (1, d("-log1p_exp(param[1])")),
        (:JAx, 1, 1) => (1, d("-log1p_exp(param[1])")),
        (:DRIFT, 2, 1) => (2, d("param[2]")), (:JAx, 2, 1) => (2, d("param[2]")),
        (:DRIFT, 1, 2) => (3, d("param[3]")), (:JAx, 1, 2) => (3, d("param[3]")),
        (:DRIFT, 2, 2) => (4, d("-log1p_exp(param[4])")),
        (:JAx, 2, 2) => (4, d("-log1p_exp(param[4])")),
        (:DIFFUSION, 1, 1) => (5, "log1p_exp(param[5])"),
        (:DIFFUSION, 2, 1) => (6, "param[6]"),
        (:DIFFUSION, 2, 2) => (7, "log1p_exp(param[7])"),
        (:CINT, 1, 1) => (8, "param[8]"), (:CINT, 2, 1) => (9, "param[9]"),
        (:LAMBDA, 2, 1) => (10, "param[10]"), (:Jy, 2, 1) => (10, "param[10]"),
        (:MANIFESTMEANS, 1, 1) => (11, "param[11]"),
        (:MANIFESTMEANS, 2, 1) => (12, "param[12]"),
        (:MANIFESTMEANS, 3, 1) => (13, "param[13]"),
        (:MANIFESTVAR, 1, 1) => (14, "log1p_exp(param[14])"),
        (:MANIFESTVAR, 2, 2) => (15, "log1p_exp(param[15])"),
        (:MANIFESTVAR, 3, 3) => (16, "log1p_exp(param[16])"),
        (:MANIFESTVAR, 3, 1) => (17, "param[17]"),
        (:T0MEANS, 1, 1) => (18, "param[18]"), (:T0MEANS, 2, 1) => (19, "param[19]"),
        (:T0VAR, 1, 1) => (20, "log1p_exp(param[20])"),
        (:T0VAR, 2, 1) => (21, "param[21]"),
        (:T0VAR, 2, 2) => (22, "log1p_exp(param[22])"),
    )
    # `perrow`: JAx[1,1] read from DRIFT at every row, as a state-dependent
    # drift's Jacobian is, which the hand route refuses.
    perrow && delete!(free, (:JAx, 1, 1))
    df = _adjoint_test_dataframe(drift=drift, jax=drift, cint=[0.1; -0.1;;],
        diffusion=[0.5 0.0; 0.1 0.4], lambda=lambda, jy=lambda,
        manifestmeans=zeros(3, 1), manifestvar=Matrix(0.3I, 3, 3),
        t0var=[1.0 0.0; 0.2 1.0], t0means=zeros(2, 1), free=free,
        predict=perrow ? Dict((:JAx, 1, 1) => "DRIFT[1,1]") : Dict())
    sp = ekf_from_data_frame(df)
    starts = Int[]
    times = Float64[]
    columns = Vector{Float64}[]
    for s in 1:6
        push!(starts, length(times) + 1)
        t = 0.0
        for r in 1:(5 + s % 4)
            push!(times, t)
            y = [0.5 * sin(1.3s + r), 0.4 * cos(0.7s - 2r), 0.3 * sin(s * r)]
            (s + r) % 5 == 0 && (y[2] = NaN)
            (s + 2r) % 7 == 0 && (y[1] = NaN)
            s == 3 && r == 2 && (y .= NaN)
            push!(columns, y)
            t += 0.3 + 0.9 * abs(sin(2.1s + 1.7r))
        end
    end
    objective = ctsem_objective(sp, starts, times, reduce(hcat, columns),
        zeros(0, length(times)), zeros(6, 0), max_timestep;
        prior_index=collect(1:22), prior_scale=fill(2.0, 22),
        subject_weights=[1.0, 2.0, 0.5, 1.0, 0.0, 1.5])
    values = [0.3, 0.15, -0.1, 0.2, -0.4, 0.2, -0.3, 0.1, -0.2, 0.6, 0.05, -0.1,
        0.15, -1.0, -0.8, -1.2, 0.1, 0.2, -0.1, 0.0, 0.3, -0.5]
    return objective, values
end

# The dual route's Hessian, whatever the default.
function _tangent_dual_hessian(objective, values)
    ContinuousTimeSEM.ctsem_set_tangent_hessian!(false)
    try
        return ctsem_hessian(objective, values)
    finally
        ContinuousTimeSEM.ctsem_set_tangent_hessian!(true)
    end
end

@testset "the hand-coded Hessian is the dual's" begin
    for max_timestep in (Inf, 0.4)
        objective, values = _tangent_test_model(max_timestep=max_timestep)
        n = length(values)
        @test ContinuousTimeSEM._tan_refusal(objective) === nothing
        H = zeros(n, n)
        # Taken by the hand route, not refused at the point.
        @test ContinuousTimeSEM._tan_hessian_columns!(H, objective, values, collect(1:n),
            nothing)
        reference = _tangent_dual_hessian(objective, values)
        @test all(isfinite, reference)
        @test (H .+ H') ./ 2 ≈ reference rtol = 1e-11
        # Unsymmetrised it is already symmetric to rounding: one exact Hessian.
        @test maximum(abs.(H .- H')) <= 1e-10 * maximum(abs.(reference))
        @test ctsem_hessian(objective, values) == (H .+ H') ./ 2
        cols = [3, 7, 13, 20]
        @test ctsem_hessian_columns(objective, values, cols) ≈ reference[:, cols] rtol = 1e-11
    end
end

@testset "an interval taken as a series goes to the dual" begin
    # A drift close to zero against the intervals: the filter integrates those
    # intervals as series, which the hand route does not cover.
    objective, values = _tangent_test_model(driftscale=1e-5)
    n = length(values)
    @test ContinuousTimeSEM._tan_refusal(objective) === nothing
    @test !ContinuousTimeSEM._tan_hessian_columns!(zeros(n, n), objective, values,
        collect(1:n), nothing)
    @test ctsem_hessian(objective, values) == _tangent_dual_hessian(objective, values)
end

@testset "what the hand route does not cover is refused" begin
    objective, values = _tangent_test_model(perrow=true)
    @test ContinuousTimeSEM._tan_refusal(objective) !== nothing
    @test ctsem_hessian(objective, values) == _tangent_dual_hessian(objective, values)
end
