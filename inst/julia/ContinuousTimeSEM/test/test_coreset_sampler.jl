# The importance sampler's coreset choice (coreset_sampler.jl) on sums it can
# compress: its training rows (`_ISRows`) filled directly, as the sampler's
# full evaluations fill them, and the choice and the weights read off them.
# The sampler's density against the full one, on a fitted model, is
# test-julia-coreset.R's.

using Test, Random, LinearAlgebra
using ContinuousTimeSEM

const CIS = ContinuousTimeSEM

# Rows into a store, as `_is_learn!` takes a batch: draws as columns of `L`.
_cis_rows(X::AbstractMatrix, y::AbstractVector) = begin
    tr = CIS._ISRows(size(X, 2))
    CIS._is_learn!(tr, (value=y, total=y, L=Matrix(transpose(X))))
    tr
end

@testset "the selection finds a small coreset where the units' terms span few dimensions" begin
    # Sixty units whose terms are each a non-negative mix of three functions
    # of the draws: their sum is a mix of the same three, so a handful of
    # units carry it, and the held-out error says so.
    rng = Random.Xoshiro(5)
    draws = 240
    basis = hcat(randn(rng, draws), randn(rng, draws) .^ 2, sin.((1:draws) ./ 9))
    mix = rand(rng, 3, 60)
    # On a log density's scale: the sum moves a few units across the draws.
    X = (basis * mix) ./ 20 .+ 1e-5 .* randn(rng, draws, 60)
    y = vec(sum(X; dims=2))
    # Chosen by its expected cost in a round of 100 draws and 30 gradients: a
    # handful of units, since three carry the sum and every unit costs.
    rd = CIS._ISRound(100, 30, 8)
    cf = CIS._is_coreset_choose(_cis_rows(X[1:200, :], y[1:200]), rd, Int[])
    @test cf.chosen
    @test count(>(0), cf.w) <= 10
    @test cf.cost < CIS._is_full_round(rd)
    held = 201:draws
    pred = [cf.c + sum(cf.w[i] * X[j, cf.S[i]] for i in eachindex(cf.S)) for j in held]
    @test maximum(abs.(pred .- y[held])) < CIS._IS_CORESET_TOL.max
    # Units that all move against the sum get zero weights, which leave
    # nothing to evaluate: never chosen, however cheap that looks, from
    # scratch or from a coreset that held them (a 200-person fit once handed
    # the engine an empty unit list).
    against = _cis_rows(50 .- X[1:200, :], y[1:200])
    @test !CIS._is_coreset_choose(against, rd, Int[]).chosen
    @test !CIS._is_coreset_choose(against, rd, [1, 2, 3]).chosen
    # Too few rows to hold any out: no choice, and why.
    few = CIS._is_coreset_choose(_cis_rows(X[1:8, :], y[1:8]), rd, Int[])
    @test !few.chosen && few.why === :held
end

@testset "a size is held to the full sum's gradient at held-out scored draws" begin
    # More parameters (30) than units (20), the draws moving in every one.
    # Units alike in how they depend on the parameters: one carries the sum's
    # gradient. Units each moving in a direction of its own: five carry none
    # of the other fifteen's.
    rng = Random.Xoshiro(11)
    N, p, D = 20, 30, 8
    Xs = randn(rng, p, D)
    alike = zeros(N, p, D)
    own = zeros(N, p, D)
    for d in 1:D
        g = randn(rng, p)
        for i in 1:N
            alike[i, :, d] .= g .* (1 + 1e-4 * randn(rng))
            own[i, i, d] = randn(rng)
        end
    end
    learnt(S) = begin
        tr = CIS._ISRows(N)
        CIS._is_learn_scores!(tr, Xs, S, zeros(D))
        tr
    end
    force(tr, units, w) = CIS._is_force_error(CIS._is_score_rows(tr, 0.0), tr,
        CIS._is_held_draws(tr.zdraw), units, w, 0.0)
    tol = CIS._IS_CORESET_TOL.force
    @test length(CIS._is_held_draws(learnt(alike).zdraw)) == 2
    @test force(learnt(alike), [1], [Float64(N)]) < tol
    @test force(learnt(own), collect(1:5), fill(4.0, 5)) > tol
    @test force(learnt(own), collect(1:N), ones(N)) < 1e-10
    # A parameter the draws do not move counts for nothing: units alike but
    # for one such parameter, where one unit is far from the sum.
    still = copy(alike)
    for i in 1:N
        still[i, 25, :] .= Float64(i)
    end
    Xstill = copy(Xs)
    Xstill[25, :] .= 0.0
    trs = CIS._ISRows(N)
    CIS._is_learn_scores!(trs, Xstill, still, zeros(D))
    @test force(trs, [1], [Float64(N)]) < tol
    trm = CIS._ISRows(N)
    CIS._is_learn_scores!(trm, Xs, still, zeros(D))
    @test force(trm, [1], [Float64(N)]) > tol
    # With no draw held out the values decide.
    tr3 = CIS._ISRows(N)
    CIS._is_learn_scores!(tr3, Xs[:, 1:3], own[:, :, 1:3], zeros(3))
    @test force(tr3, [1], [1.0]) == 0.0
end

@testset "score rows fit the weights that values alone leave open" begin
    # Twelve units whose terms are each a non-negative mix of three functions
    # of a two-parameter draw, and their scores. Four draws' values cannot fit
    # six weights; their scores add two rows each, and with them the coreset
    # that reproduces the sum is found and predicts held-out draws exactly.
    # The rows are in units of each parameter's spread, so rescaling a
    # parameter changes nothing.
    rng = Random.Xoshiro(7)
    f(th) = hcat(th[:, 1], th[:, 2] .^ 2, sin.(th[:, 1]) .* th[:, 2])
    df(t) = ([1.0, 0.0, cos(t[1]) * t[2]], [0.0, 2 * t[2], sin(t[1])])
    mix = zeros(3, 12)
    mix[:, 1:6] .= rand(rng, 3, 6)
    ytrue = vec(sum(mix; dims=2))
    th = randn(rng, 4, 2)
    X = f(th) * mix
    y = f(th) * ytrue
    stored(scale) = begin
        tr = _cis_rows(X, y)
        for r in 1:4
            S = zeros(12, 2, 1)
            for j in 1:2
                S[:, j, 1] .= vec(transpose(df(th[r, :])[j]) * mix) ./ scale[j]
            end
            CIS._is_learn_scores!(tr, reshape(th[r, :] .* scale, 2, 1), S, [maximum(y)])
        end
        tr
    end
    w, c = CIS._is_weights(stored([1.0, 1.0]), collect(1:6))
    heldth = randn(rng, 10, 2)
    @test f(heldth) * mix[:, 1:6] * w .+ c ≈ f(heldth) * ytrue rtol = 1e-6
    # The second parameter on a scale a thousand times larger: the same weights.
    wbig, _ = CIS._is_weights(stored([1.0, 1000.0]), collect(1:6))
    @test wbig ≈ w rtol = 1e-6
end

@testset "a check reads the coreset's top half and an even spread of the rest" begin
    v = collect(1.0:20.0)
    pick = CIS._is_checkpick(v)
    @test length(pick) == 8
    @test sort(pick[1:4]) == [17, 18, 19, 20]
    @test all(p -> 1 <= p <= 16, pick[5:8])
    @test CIS._is_checkpick([3.0, 1.0]) == [1, 2]
    @test CIS._is_checkpick([5.0]) == [1]
end

@testset "relevance is measured from the run's best, or the newest batch's" begin
    L = ones(2, 6)
    learn!(tr, y) = CIS._is_learn!(tr, (value=y, total=y, L=L[:, 1:length(y)]))
    for batchtop in (false, true)
        tr = CIS._ISRows(2; batchtop=batchtop)
        CIS._is_batch_start!(tr)
        learn!(tr, [0.0, -1.0, -2.0])
        @test CIS._is_top(tr) == 0.0
        CIS._is_batch_start!(tr)
        learn!(tr, [-100.0, -90.0])
        @test CIS._is_top(tr) == (batchtop ? -90.0 : 0.0)
        # A batch that has learned nothing yet reads every row.
        CIS._is_batch_start!(tr)
        @test CIS._is_top(tr) == 0.0
    end
end
