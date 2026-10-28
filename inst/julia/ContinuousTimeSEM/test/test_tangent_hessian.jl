using DataFrames, ForwardDiff, LinearAlgebra

# The hand-coded second order (`tangent_hessian.jl`) against the Hessian's dual
# route, its oracle: the two compute the same derivative of the same filter, so
# they agree to rounding. One model exercises what the hand route covers --
# every drift cell, a cross-loading, correlated diffusion and measurement
# errors, missing values and a row with none observed, substeps, unequal and
# zero subject weights, a prior, two chunks of columns, TI predictors whose
# effects give each subject its own parameter layer, a TD predictor, stationary
# initial moments, the same cells in discrete time, and cells the filter
# evaluates by row, reading the state.

isdefined(@__MODULE__, :_adjoint_test_dataframe) ||
    include(joinpath(@__DIR__, "adjoint_fixtures.jl"))

function _tangent_test_model(; driftscale=1.0, max_timestep=Inf, perrow=false, ti=false,
    td=false, stationary=false, continuous=true, rows=false, relin=false, plan=-1,
    moment=false, cat=false, timiss=false, elements=false)
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
    # `rows`: cells the filter evaluates by row, reading the state -- a bilinear
    # drift, its Jacobian, a diffusion and a manifest mean -- through free PARS
    # cells (parameters 25 and 26), beside the TD predictor.
    predict = perrow ? Dict((:JAx, 1, 1) => "DRIFT[1,1]") : Dict{Tuple{Symbol,Int,Int},String}()
    update = Dict{Tuple{Symbol,Int,Int},String}()
    if rows
        td = true
        for key in ((:JAx, 1, 1), (:DRIFT, 1, 2), (:JAx, 1, 2), (:DIFFUSION, 2, 2),
            (:MANIFESTMEANS, 3, 1))
            delete!(free, key)
        end
        free[(:PARS, 1, 1)] = (25, "param[25]")
        free[(:PARS, 2, 1)] = (26, "param[26]")
        predict[(:DRIFT, 1, 2)] = "PARS[1,1] * state[1]"
        predict[(:JAx, 1, 1)] = "DRIFT[1,1] + PARS[1,1] * state[2]"
        predict[(:JAx, 1, 2)] = "PARS[1,1] * state[1]"
        predict[(:DIFFUSION, 2, 2)] = "log1p_exp(PARS[2,1]) * (1 + 0.1 * state[1]^2)"
        update[(:MANIFESTMEANS, 3, 1)] = "PARS[2,1] + 0.2 * state[1]^2"
    end
    # `relin`: a loading and its Jacobian reading the state, so the filter
    # iterates each row's update to its mode; `plan`, a frozen count of steps
    # a row (-1 to the mode).
    # `moment`: the first drift cell cubic in its state, its Jacobian bending, so
    # the filter moment-matches each prediction.
    if moment
        delete!(free, (:DRIFT, 1, 1))
        predict[(:DRIFT, 1, 1)] = "-log1p_exp(PARS[2,1]) * (1 + 0.2 * state[1]^2)"
        predict[(:JAx, 1, 1)] =
            "-log1p_exp(PARS[2,1]) * (1 + 0.6 * state[1]^2) + PARS[1,1] * state[2]"
    end
    if relin
        delete!(free, (:LAMBDA, 2, 1))
        delete!(free, (:Jy, 2, 1))
        update[(:LAMBDA, 2, 1)] = "PARS[1,1] + 0.3 * state[1]"
        update[(:Jy, 2, 1)] = "PARS[1,1] + 0.6 * state[1]"
    end
    # `cat`: the second manifest ordinal in three categories, its thresholds
    # parameters 27 and 28, and the third binary.
    if cat
        free[(:THRESHOLDS, 2, 1)] = (27, "param[27]")
        free[(:THRESHOLDS, 2, 2)] = (28, "log1p_exp(param[28])")
    end
    # `td`: a TD predictor whose effects on both states are parameters 23 and 24.
    if td
        free[(:TDPREDEFFECT, 1, 1)] = (23, "param[23]")
        free[(:TDPREDEFFECT, 2, 1)] = (24, "param[24]")
    end
    df = _adjoint_test_dataframe(drift=drift, jax=drift, cint=[0.1; -0.1;;],
        diffusion=[0.5 0.0; 0.1 0.4], lambda=lambda, jy=lambda,
        manifestmeans=zeros(3, 1), manifestvar=Matrix(0.3I, 3, 3),
        t0var=[1.0 0.0; 0.2 1.0], t0means=zeros(2, 1), free=free,
        tdpredeffect=td ? zeros(2, 1) : nothing, jtd=td ? [1.0 0.0; 0.0 1.0] : nothing,
        pars=rows ? zeros(2, 1) : zeros(1, 1), predict=predict, update=update,
        extra=cat ? [:THRESHOLDS => zeros(3, 2)] : Pair{Symbol,Matrix{Float64}}[])
    # `ti`: two TI predictors, with effects on a drift cell, a diffusion sd and
    # correlation, an intercept, a loading, a mean, a measurement sd, an
    # initial mean and an initial sd -- each subject a layer of its own.
    tieffects = DataFrame(parameter=[1, 5, 6, 8, 10, 12, 14, 18, 20],
        predictor=[1, 2, 1, 2, 1, 2, 1, 2, 1], coefficient=collect(23:31))
    # `stationary`: the first row's prior the stationary moments; `continuous`
    # false, the same cells as a discrete-time model.
    sp = ekf_from_data_frame(df, ti ? tieffects : DataFrame(parameter=Int[],
        predictor=Int[], coefficient=Int[]); continuous_time=continuous,
        stationary=stationary, manifesttype=cat ? [0, 2, 1] : Int[],
        ncategories=cat ? [0, 3, 2] : Int[],
        # `elements`: the first latent's initial mean and variance stationary,
        # over the closure the second drives into it.
        stationary_mean=elements ? [1, 0] : Int[], stationary_block=elements ? [1, 0] : Int[],
        stationary_closure=elements ? [1, 1] : Int[])
    npar = ti ? (timiss ? 33 : 31) : cat ? 28 : rows ? 26 : td ? 24 : 22
    starts = Int[]
    times = Float64[]
    columns = Vector{Float64}[]
    tdcolumns = Vector{Float64}[]
    for s in 1:6
        push!(starts, length(times) + 1)
        t = 0.0
        for r in 1:(5 + s % 4)
            push!(times, t)
            y = [0.5 * sin(1.3s + r), 0.4 * cos(0.7s - 2r), 0.3 * sin(s * r)]
            if cat
                y[2] = 1.0 + (y[2] > -0.15) + (y[2] > 0.2)
                y[3] = y[3] > 0 ? 1.0 : 0.0
            end
            (s + r) % 5 == 0 && (y[2] = NaN)
            (s + 2r) % 7 == 0 && (y[1] = NaN)
            s == 3 && r == 2 && (y .= NaN)
            push!(columns, y)
            push!(tdcolumns, [r % 3 == 1 ? 1.0 : 0.4 * sin(s + r)])
            t += 0.3 + 0.9 * abs(sin(2.1s + 1.7r))
        end
    end
    tipreds = ti ? [0.3 -1.1; -0.8 0.4; 1.2 0.9; 0.0 -0.5; 0.5 1.4; -1.3 -0.2] :
        zeros(6, 0)
    # `timiss`: the first predictor missing for subjects 2 and 5, its values
    # sampled as parameters 32 and 33 (`TIMissingRecipe`).
    timiss && (tipreds[2, 1] = 99999.0; tipreds[5, 1] = 99999.0)
    missing_kw = timiss ? (ti_missing_subject=[2, 5], ti_missing_predictor=[1, 1],
        ti_missing_parameter=[32, 33], ti_missing_mu=[0.1, 0.1],
        ti_missing_sigma=[0.8, 0.8]) : (;)
    objective = ctsem_objective(sp, starts, times, reduce(hcat, columns),
        td ? reduce(hcat, tdcolumns) : zeros(0, length(times)), tipreds, max_timestep;
        prior_index=collect(1:npar), prior_scale=fill(2.0, npar),
        subject_weights=[1.0, 2.0, 0.5, 1.0, 0.0, 1.5], missing_kw...)
    plan >= 0 && foreach(so -> append!(so.relin, fill(plan, length(so.timesteps))),
        objective.subject_objectives)
    values = [0.3, 0.15, -0.1, 0.2, -0.4, 0.2, -0.3, 0.1, -0.2, 0.6, 0.05, -0.1,
        0.15, -1.0, -0.8, -1.2, 0.1, 0.2, -0.1, 0.0, 0.3, -0.5]
    ti && append!(values, [0.1, -0.15, 0.2, 0.05, -0.1, 0.12, 0.08, -0.2, 0.15])
    timiss && append!(values, [0.4, -0.6])
    td && append!(values, [0.7, -0.4])
    rows && append!(values, [0.3, -0.2])
    cat && (append!(values, zeros(28 - length(values))); values[27:28] = [-0.3, 0.2])
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
    for opts in ((;), (max_timestep=0.4,), (ti=true,), (max_timestep=0.4, td=true),
        (stationary=true,), (continuous=false, td=true), (rows=true,),
        (rows=true, max_timestep=0.4), (rows=true, continuous=false), (perrow=true,),
        (rows=true, relin=true), (rows=true, relin=true, plan=2),
        (rows=true, moment=true), (rows=true, moment=true, relin=true, max_timestep=0.4),
        (rows=true, moment=true, continuous=false), (cat=true,), (cat=true, td=true),
        (cat=true, rows=true, relin=true), (cat=true, continuous=false, max_timestep=0.4),
        (stationary=true, perrow=true), (elements=true,), (elements=true, perrow=true),
        (ti=true, timiss=true))
        objective, values = _tangent_test_model(; opts...)
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

# The augmented route's shape: two latents and a third state carrying a
# random intercept on the first's drift, static, with no diffusion and its
# population variance a block of its own (`RAWPOPVAR`, `_place_population_block!`);
# the intercept's solve then runs over the two latents and the noise's over
# the diffusing ones (`affine_dim`, `diffusion_state_indices`).
function _tangent_aug_model(; continuous=true, max_timestep=Inf, driftscale=1.0,
    stationary=false)
    carry = continuous ? 0.0 : 1.0
    drift = [-0.6 0.2 1.0; 0.1 -0.5 0.0; 0.0 0.0 carry]
    lambda = [1.0 0.0 0.0; 0.7 0.0 0.0; 0.0 1.0 0.0]
    d(e) = driftscale == 1 ? e : "$(driftscale) * ($(e))"
    free = Dict(
        (:DRIFT, 1, 1) => (1, d("-log1p_exp(param[1])")),
        (:JAx, 1, 1) => (1, d("-log1p_exp(param[1])")),
        (:DRIFT, 2, 1) => (2, d("param[2]")), (:JAx, 2, 1) => (2, d("param[2]")),
        (:DRIFT, 2, 2) => (3, d("-log1p_exp(param[3])")),
        (:JAx, 2, 2) => (3, d("-log1p_exp(param[3])")),
        (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
        (:DIFFUSION, 2, 1) => (5, "param[5]"),
        (:DIFFUSION, 2, 2) => (6, "log1p_exp(param[6])"),
        (:CINT, 2, 1) => (7, "param[7]"),
        (:LAMBDA, 2, 1) => (8, "param[8]"), (:Jy, 2, 1) => (8, "param[8]"),
        (:MANIFESTMEANS, 1, 1) => (9, "param[9]"),
        (:MANIFESTMEANS, 3, 1) => (10, "param[10]"),
        (:MANIFESTVAR, 1, 1) => (11, "log1p_exp(param[11])"),
        (:MANIFESTVAR, 2, 2) => (12, "log1p_exp(param[12])"),
        (:MANIFESTVAR, 3, 3) => (13, "log1p_exp(param[13])"),
        (:T0MEANS, 1, 1) => (14, "param[14]"),
        (:T0VAR, 1, 1) => (15, "log1p_exp(param[15])"),
        (:T0VAR, 2, 2) => (16, "log1p_exp(param[16])"),
        (:RAWPOPVAR, 1, 1) => (17, "log1p_exp(param[17])"),
        (:T0MEANS, 3, 1) => (18, "param[18]"),
    )
    df = _adjoint_test_dataframe(drift=drift, jax=drift, cint=zeros(3, 1),
        diffusion=zeros(3, 3), lambda=lambda, jy=lambda, manifestmeans=zeros(3, 1),
        manifestvar=Matrix(0.3I, 3, 3), t0var=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0],
        t0means=zeros(3, 1), free=free, extra=[:RAWPOPVAR => zeros(1, 1)])
    sp = ekf_from_data_frame(df, DataFrame(parameter=Int[], predictor=Int[],
        coefficient=Int[]), [1, 2]; continuous_time=continuous, affine_dim=2,
        population_indices=[3], population_scale=[1.0], stationary=stationary)
    starts = Int[]
    times = Float64[]
    columns = Vector{Float64}[]
    for s in 1:6
        push!(starts, length(times) + 1)
        t = 0.0
        for r in 1:(5 + s % 4)
            push!(times, t)
            y = [0.5 * sin(1.3s + r) + 0.3 * sin(s), 0.4 * cos(0.7s - 2r), 0.3 * sin(s * r)]
            (s + r) % 5 == 0 && (y[2] = NaN)
            push!(columns, y)
            t += 0.3 + 0.9 * abs(sin(2.1s + 1.7r))
        end
    end
    objective = ctsem_objective(sp, starts, times, reduce(hcat, columns),
        zeros(0, length(times)), zeros(6, 0), max_timestep;
        prior_index=collect(1:18), prior_scale=fill(2.0, 18),
        subject_weights=[1.0, 2.0, 0.5, 1.0, 0.0, 1.5])
    values = [0.3, 0.15, 0.2, -0.4, 0.2, 0.1, -0.2, 0.6, 0.05, -0.1, -1.0, -0.8, -1.2,
        0.1, 0.0, -0.5, -0.3, 0.2]
    return objective, values
end

@testset "the augmented route's static carriers and population block" begin
    for opts in ((;), (max_timestep=0.4,), (continuous=false,), (driftscale=0.02,),
        (stationary=true,))
        objective, values = _tangent_aug_model(; opts...)
        n = length(values)
        @test ContinuousTimeSEM._tan_refusal(objective) === nothing
        H = zeros(n, n)
        @test ContinuousTimeSEM._tan_hessian_columns!(H, objective, values, collect(1:n),
            nothing)
        reference = _tangent_dual_hessian(objective, values)
        @test all(isfinite, reference)
        @test (H .+ H') ./ 2 ≈ reference rtol = 1e-10
    end
end

@testset "the row scores by the hand-coded forward are the filter's" begin
    # `ctsem_row_scores` along 15 directions, two chunks of lanes, against the
    # filter at the Hessian's dual, its oracle.
    for (make, opts) in ((_tangent_test_model, (;)), (_tangent_test_model, (ti=true,)),
        (_tangent_test_model, (rows=true, relin=true, max_timestep=0.4)),
        (_tangent_test_model, (cat=true, moment=true, rows=true)),
        (_tangent_aug_model, (;)))
        objective, values = make(; opts...)
        dirs = [sin(1.7i + 0.9j) for i in eachindex(values), j in 1:15]
        hand = ctsem_row_scores(objective, values, dirs)
        rows = [size(s.data, 2) for s in objective.subject_objectives]
        H = fill(NaN, sum(rows), 15)
        @test ContinuousTimeSEM._tan_row_scores!(H, objective, collect(values), dirs, rows)
        ContinuousTimeSEM.ctsem_set_tangent_hessian!(false)
        dual = try
            ctsem_row_scores(objective, values, dirs)
        finally
            ContinuousTimeSEM.ctsem_set_tangent_hessian!(true)
        end
        @test all(isfinite, dual)
        @test hand ≈ dual rtol = 1e-10
    end
end

@testset "intervals taken as series" begin
    # A drift close to zero against the intervals: the filter integrates every
    # interval's intercept and noise as series. At 0.02, with substeps, some
    # intervals take them and some the closed forms.
    for opts in ((driftscale=1e-5,), (driftscale=0.02, max_timestep=0.4))
        objective, values = _tangent_test_model(; opts...)
        n = length(values)
        counts = ContinuousTimeSEM._CTSEM_OPCOUNT
        before = (counts.series_intercept[], counts.series_noise[])
        reference = _tangent_dual_hessian(objective, values)
        # The filter took both series somewhere.
        @test counts.series_intercept[] > before[1]
        @test counts.series_noise[] > before[2]
        H = zeros(n, n)
        @test ContinuousTimeSEM._tan_hessian_columns!(H, objective, values, collect(1:n),
            nothing)
        @test (H .+ H') ./ 2 ≈ reference rtol = 1e-10
    end
end

@testset "the test models take the filter's routes they are meant to" begin
    # Cells by row without iterating the update; with a loading's Jacobian
    # reading the state, iterated; neither moment-matched.
    for (relin, moment) in ((false, false), (true, false), (false, true))
        objective, _ = _tangent_test_model(rows=true, relin=relin, moment=moment)
        ws = ContinuousTimeSEM._init_continuous_ekf_workspace(Float64, objective.params)
        @test ws.relinearise == relin
        @test ws.moment == moment
    end
end
