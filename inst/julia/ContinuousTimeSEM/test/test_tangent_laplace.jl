# The hand-coded Laplace unit curvature (`tangent_laplace.jl`) against the
# member gradient's Jacobian at the curvature's dual, its oracle
# (`ctsem_set_tangent_hessian!(false)`): the dense and blocked routes, one to
# three levels, effects on the means and on the dynamics and noise, and a
# drift evaluated by row.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

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
