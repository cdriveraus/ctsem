# The hand-coded Laplace unit curvature (`tangent_laplace.jl`) against the
# member gradient's Jacobian at the curvature's dual, its oracle
# (`ctsem_set_tangent_hessian!(false)`): the dense and blocked routes, one to
# three levels, effects on the means and on the dynamics and noise, both at
# once by the split route, a curved mean, and a drift evaluated by row.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

# The linear model with its initial mean's transform curved: an effect there
# moves the innovations nonlinearly, so Gauss-Newton is not its curvature, and
# the objective says so as R's would, `affine` naming only the intercept.
const _TAN_CURVED_MEAN_OBJECTIVE = let
    df = _laplace_test_dataframe(
        drift=[-0.5;;], jax=[-0.5;;], cint=[0.0;;], diffusion=[0.2;;],
        lambda=[1.0;;], jy=[1.0;;], manifestmeans=[0.0;;], manifestvar=[0.3;;],
        t0var=[0.4;;], t0means=[0.0;;],
        free=Dict(
            (:T0MEANS, 1, 1) => (1, "2 * log1p_exp(param[1])"),
            (:CINT, 1, 1) => (2, "param[2]"),
            (:DRIFT, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:JAx, 1, 1) => (3, "-log1p_exp(param[3])"),
            (:DIFFUSION, 1, 1) => (4, "log1p_exp(param[4])"),
            (:MANIFESTMEANS, 1, 1) => (5, "param[5]"),
        ),
    )
    subject_starts, times, data = _laplace_test_data(6, 1)
    ctsem_objective(ekf_from_data_frame(df), subject_starts, times, data)
end

@testset "hand-coded Laplace curvature against the dual" begin
    C = ContinuousTimeSEM
    function curvatures(laplace, values, seed)
        Ls = C._laplace_popchols(values, laplace.spec)
        rng = Random.MersenneTwister(seed)
        out = Tuple{Matrix{Float64},Matrix{Float64}}[]
        for U in eachindex(laplace.units.members)
            d = length(laplace.units.blocks[U]) == 0 ? 0 :
                sum(b -> b.size, laplace.units.blocks[U])
            d == 0 && continue
            u = 0.4 .* randn(rng, d)
            dense(on) = (C.ctsem_set_tangent_hessian!(on);
                C._laplace_block_dense(C._laplace_unit_curvature(laplace, U, values, Ls, u),
                    laplace.units.blocks[U], d))
            push!(out, (dense(true), dense(false)))
        end
        C.ctsem_set_tangent_hessian!(true)
        return out
    end
    function served(laplace, values)
        ctx = get(laplace.workspaces[C._laplace_slot()],
            (:tan_curvature, length(values), C._HYBRID_PROFILE), nothing)
        return ctx !== nothing
    end
    cases = [
        ("means", () -> _fresh_linear()),
        # Every effect in the means affinely: Gauss-Newton for the whole member.
        ("means, quadratic", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [1, 2], [6, 7], [8], [1.0, 1.0]; quadratic=true), copy(_LAPLACE_LINEAR_VALUES))),
        # Two effects on the dynamics and noise beside three in the means: the
        # split, one full sweep and one forward (`_tan_laplace_split_member!`).
        ("means beside dynamics, split", () -> (ctsem_laplace_objective(
            _LAPLACE_LINEAR_OBJECTIVE, [1, 2, 5, 3, 4], collect(6:10), collect(11:20),
            ones(5); affine_level=[1, 1, 1], affine_effect=[1, 2, 3]),
            vcat(_LAPLACE_LINEAR_VALUES[1:5], [-0.3, -0.15, -0.2, -0.4, -0.35],
                0.1 .* sin.(1:10)))),
        # An intercept beside an initial mean whose transform is curved: the
        # objective's `affine` keeps the member off Gauss-Newton.
        ("curved mean", () -> (ctsem_laplace_objective(_TAN_CURVED_MEAN_OBJECTIVE,
            [1, 2], [6, 7], [8], [1.0, 1.0]; affine_level=[1], affine_effect=[2]),
            copy(_LAPLACE_LINEAR_VALUES))),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("drift by row", () -> _fresh_nonlinear()),
        ("two levels", () -> _fresh_twolevel()),
        ("three levels", () -> _fresh_threelevel()),
    ]
    threshold = C._LAPLACE_BLOCK_THRESHOLD[]
    try
        for (name, make) in cases, blocked in (false, true)
            C.ctsem_set_block_threshold!(blocked ? 1 : threshold)
            laplace, values = make()
            pairs = curvatures(laplace, values, 3)
            @testset "$name, $(blocked ? "blocked" : "dense")" begin
                @test served(laplace, values)
                @test !isempty(pairs)
                for (hand, dual) in pairs
                    @test all(isfinite, hand)
                    @test maximum(abs, hand .- dual) <= 1e-9 * max(1.0, maximum(abs, dual))
                end
            end
        end
    finally
        C.ctsem_set_block_threshold!(threshold)
        C.ctsem_set_tangent_hessian!(true)
    end
end

@testset "hand-coded seeded sweeps against the dual" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("drift by row", () -> _fresh_nonlinear()),
        ("three levels", () -> _fresh_threelevel()),
    ]
    try
        for (name, make) in cases
            laplace, values = make()
            Ls = C._laplace_popchols(values, laplace.spec)
            rng = Random.MersenneTwister(7)
            @testset "$name" begin
                for U in eachindex(laplace.units.members)
                    d = sum(b -> b.size, laplace.units.blocks[U]; init=0)
                    u = 0.4 .* randn(rng, d)
                    nm = length(laplace.units.members[U])
                    for d1 in (randn(rng, length(values)), zeros(length(values)))
                        C.ctsem_set_tangent_hessian!(true)
                        hand = C._laplace_unit_seeded_gradient(laplace, U, values, Ls, u,
                            1:nm, d1)
                        h0, h1 = copy(hand.d0), copy(hand.d1c)
                        C.ctsem_set_tangent_hessian!(false)
                        dual = C._laplace_unit_seeded_gradient(laplace, U, values, Ls, u,
                            1:nm, d1)
                        @test hand.ok && dual.ok
                        @test maximum(abs, h0 .- dual.d0) <= 1e-9 * max(1.0, maximum(abs, dual.d0))
                        @test maximum(abs, h1 .- dual.d1c) <= 1e-9 * max(1.0, maximum(abs, dual.d1c))
                    end
                end
            end
        end
    finally
        C.ctsem_set_tangent_hessian!(true)
    end
end

@testset "hand-coded sequential predictions against the dual" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("drift by row", () -> _fresh_nonlinear()),
    ]
    try
        for (name, make) in cases
            laplace, values = make()
            C.ctsem_set_tangent_hessian!(true)
            hand = C.ctsem_laplace_sequential(laplace, values)
            served = get(laplace.workspaces[C._laplace_slot()],
                (:tan_curvature, length(values), C._HYBRID_OFF), nothing) !== nothing
            C.ctsem_set_tangent_hessian!(false)
            dual = C.ctsem_laplace_sequential(laplace, values)
            @testset "$name" begin
                @test served
                for f in (:eta, :etacov, :y, :ycov, :llrow, :effect, :effectcov)
                    a, b = getproperty(hand, f), getproperty(dual, f)
                    @test isequal(isnan.(a), isnan.(b))
                    ok = .!isnan.(b)
                    @test any(ok)
                    @test maximum(abs, a[ok] .- b[ok]) <= 1e-9 * max(1.0, maximum(abs, b[ok]))
                end
            end
        end
    finally
        C.ctsem_set_tangent_hessian!(true)
    end
end

@testset "hand-coded mode Jacobian and gated adjoint against the dual" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("drift by row", () -> _fresh_nonlinear()),
        ("two levels", () -> _fresh_twolevel()),
    ]
    rel(a, b) = maximum(abs, a .- b) / max(1.0, maximum(abs, b))
    try
        for (name, make) in cases
            laplace, values = make()
            C.ctsem_laplace_evaluate(laplace, values; gradient=false)
            run(on) = (C.ctsem_set_tangent_hessian!(on);
                C.ctsem_laplace_mode_jacobian(laplace, values))
            hand, dual = run(true), run(false)
            adj(on) = (C.ctsem_set_tangent_hessian!(on);
                [C._laplace_gated_unit_adjoint(laplace, U, values)
                 for U in eachindex(laplace.modes)])
            ahand, adual = adj(true), adj(false)
            @testset "$name" begin
                @test length(hand) == length(dual)
                for (a, b) in zip(hand, dual)
                    @test all(isfinite, a)
                    @test rel(a, b) <= 1e-9
                end
                @test any(!isnothing, adual)
                for (a, b) in zip(ahand, adual)
                    @test (a === nothing) == (b === nothing)
                    b === nothing || @test rel(a, b) <= 1e-8
                end
            end
        end
    finally
        C.ctsem_set_tangent_hessian!(true)
    end
end

@testset "hand-coded trace jets against the jet" begin
    C = ContinuousTimeSEM
    cases = [
        ("means", () -> _fresh_linear()),
        ("dynamics and noise", () -> (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE,
            [3, 4], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("means beside dynamics", () -> (ctsem_laplace_objective(
            _LAPLACE_LINEAR_OBJECTIVE, [1, 2, 5, 3, 4], collect(6:10), collect(11:20),
            ones(5)), vcat(_LAPLACE_LINEAR_VALUES[1:5], [-0.3, -0.15, -0.2, -0.4, -0.35],
            0.1 .* sin.(1:10)))),
        ("curved mean", () -> (ctsem_laplace_objective(_TAN_CURVED_MEAN_OBJECTIVE,
            [1, 2], [6, 7], [8], [1.0, 1.0]), copy(_LAPLACE_LINEAR_VALUES))),
        ("two levels", () -> _fresh_twolevel()),
        ("three levels", () -> _fresh_threelevel()),
        # Cells by row: the engine's own code at the jet, a chunk at a time.
        ("drift by row", () -> _fresh_nonlinear()),
    ]
    try
        for (name, make) in cases
            laplace, values = make()
            Ls = C._laplace_popchols(values, laplace.spec)
            rng = Random.MersenneTwister(11)
            npar = length(values)
            @testset "$name" begin
                for U in eachindex(laplace.units.members)
                    d = sum(b -> b.size, laplace.units.blocks[U]; init=0)
                    u = 0.4 .* randn(rng, d)
                    nm = length(laplace.units.members[U])
                    # One direction (the narrow jet), three, and seven (two
                    # chunks of the sweep's width).
                    for nd in (1, 3, 7)
                        D = zeros(npar, C._CTSEM_JET_WIDTH)
                        D[:, 1:nd] .= randn(rng, npar, nd)
                        C.ctsem_set_tangent_hessian!(true)
                        hand = C._laplace_unit_jet_gradient(laplace, U, values, Ls, u,
                            1:nm, D, nd)
                        h0, hS = copy(hand.d0), copy(hand.dS)
                        C.ctsem_set_tangent_hessian!(false)
                        jet = C._laplace_unit_jet_gradient(laplace, U, values, Ls, u,
                            1:nm, D, nd)
                        @test hand.ok && jet.ok
                        @test maximum(abs, h0 .- jet.d0) <=
                            1e-9 * max(1.0, maximum(abs, jet.d0))
                        @test maximum(abs, hS .- jet.dS) <=
                            1e-9 * max(1.0, maximum(abs, jet.dS))
                    end
                end
                served = get(laplace.workspaces[C._laplace_slot()],
                    (:tan_jet, npar, C.CTSEMTraceJet{C._CTSEM_JET_WIDTH}, C._HYBRID_PROFILE),
                    nothing) !== nothing
                @test served
            end
        end
    finally
        C.ctsem_set_tangent_hessian!(true)
    end
end
