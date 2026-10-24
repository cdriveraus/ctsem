using DataFrames, LinearAlgebra, Random, Statistics

# SAEM (saem.jl). Every check is against a closed form, not against Laplace's
# answer on a model where Laplace might share a mistake: on the linear
# fixtures the random effects enter through identity transforms of a
# linear-Gaussian model, so a unit's integrand is exactly Gaussian in its
# latent vector. There
#   * the conditional of u given the data is exactly N(mode, M^-1), so the
#     E-step's draws must have that mean and covariance;
#   * the Fisher identity's expectation of the complete-data score equals the
#     marginal gradient, which Laplace computes exactly;
#   * SAEM's fixed point is the marginal posterior mode, which the Laplace
#     optimiser finds exactly.
# Moments are checked at two draw counts, because a wrong sampler returns a
# plausible wrong answer whose error does not shrink.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

const _S = ContinuousTimeSEM

# The linear objective with a N(0, 1) prior on every raw parameter, sharing its
# compiled subject objectives, so a posterior mode is interior on six subjects.
function _saem_prior_objective(npar)
    o = _LAPLACE_LINEAR_OBJECTIVE
    _S.CTSEMObjective(o.params, o.subject_objectives, nothing, collect(1:npar),
        ones(npar), 1.0, Int[], Float64[], Float64[], Int[])
end

"""E-step only, at fixed theta: `n` sweeps of unit `U`, the draws as columns."""
function _saem_estep_draws(laplace, values, U, n; seed=1, independence=false)
    st = ctsem_saem_init(laplace, values; seed=seed)
    Ls = _S._laplace_popchols(st.theta, laplace.spec)
    d = length(st.u[U][1])
    draws = zeros(d, n)
    for k in 1:n
        st.iteration += 1
        _S._laplace_parallel(laplace, [U]) do U
            (k % 25 == 0) && _S._saem_refresh!(st, laplace, U, Ls)
            _S._saem_sweep!(st, laplace, U, 1, Ls, 1, 1 / (1 + k)^0.6;
                independence=independence)
            true
        end
        draws[:, k] = st.u[U][1]
    end
    return draws, st
end

cov2cor(C) = C ./ (sqrt.(diag(C)) * sqrt.(diag(C))')

"""The exact conditional of unit `U`'s latent vector: its mode and M^-1."""
function _saem_exact_conditional(laplace, values, U)
    ctsem_laplace_evaluate(laplace, values; gradient=false)
    theta = collect(Float64, values)
    Ls = _S._laplace_popchols(theta, laplace.spec)
    M = nothing
    _S._laplace_parallel(laplace, [U]) do U
        M = _S._laplace_unit_curvature(laplace, U, theta, Ls, laplace.modes[U])
        true
    end
    d = laplace.units.dims[U]
    dense = _S._laplace_block_dense(M, laplace.units.blocks[U], d)
    return copy(laplace.modes[U]), inv(Symmetric(dense))
end

@testset "SAEM: member log likelihoods are the inner objective's" begin
    for fresh in (_fresh_linear, _fresh_twolevel, _fresh_threelevel)
        laplace, values = fresh()
        st = ctsem_saem_init(laplace, values; seed=3)
        for U in eachindex(laplace.units.members)
            for c in 1:st.chains
                @test sum(st.ll[U][c]) ≈ _unit_loglik(laplace, U, values, st.u[U][c]) rtol = 1e-12
            end
        end
    end
end

@testset "SAEM: E-step draws the exact conditional ($name, $kernel)" for (name, fresh) in (
        ("one level", _fresh_linear), ("two levels", _fresh_twolevel),
        ("three levels", _fresh_threelevel)), kernel in ("random walk", "with independence moves")
    laplace, values = fresh()
    U = 1
    mode, Sigma = _saem_exact_conditional(laplace, values, U)
    sd = sqrt.(diag(Sigma))
    errs = map((4000, 16000)) do n
        draws, st = _saem_estep_draws(laplace, values, U, n + 500; seed=7,
            independence=(kernel != "random walk"))
        x = draws[:, 501:end]
        m = vec(mean(x; dims=2))
        C = cov(x; dims=2)
        (mean=maximum(abs.(m .- mode) ./ sd),
         var=maximum(abs.(diag(C) ./ diag(Sigma) .- 1)),
         cor=maximum(abs.(cov2cor(C) .- cov2cor(Sigma))),
         acc=_S._saem_acceptance(st))
    end
    @info "SAEM E-step moments ($name, $kernel)" errs
    # At the larger count: means within a tenth of a posterior sd, variances
    # within 15%, correlations within 0.1 -- and none of them worse than at
    # the smaller count by more than its own sampling error would allow.
    @test errs[2].mean < 0.1
    @test errs[2].var < 0.15
    @test errs[2].cor < 0.1
    @test errs[2].mean < errs[1].mean + 0.05
    @test 0.1 < errs[2].acc < 0.7
end


@testset "SAEM: the Fisher-identity gradient is the marginal gradient" begin
    laplace, values = _fresh_twolevel()
    exact = ctsem_laplace_evaluate(laplace, values; gradient=true).gradient
    st = ctsem_saem_init(laplace, values; seed=11)
    theta = collect(Float64, values)
    Ls = _S._laplace_popchols(theta, laplace.spec)
    dL = _S._laplace_level_chol_derivatives(theta, laplace.spec)
    positions = [_S._laplace_level_positions(laplace.spec, l)
                 for l in eachindex(laplace.spec.levels)]
    conditionals = [_saem_exact_conditional(laplace, values, U)
                    for U in eachindex(laplace.units.members)]
    rng = Random.Xoshiro(5)
    K = 3000
    G = zeros(length(theta), K)
    for k in 1:K
        for U in eachindex(laplace.units.members)
            mode, Sigma = conditionals[U]
            st.u[U][1] = mode .+ cholesky(Symmetric(Sigma)).L * randn(rng, length(mode))
            _S._laplace_parallel(laplace, [U]) do U
                _S._saem_unit_scores!(st, laplace, U, 1, Ls, dL, positions)
            end
            G[:, k] .+= vec(sum(st.scores[U][1]; dims=2))
        end
    end
    g = vec(mean(G; dims=2))
    se = vec(std(G; dims=2)) ./ sqrt(K)
    @test all(abs.(g .- exact) .<= 4 .* se .+ 1e-8)
end

@testset "SAEM: fixed point is the marginal posterior mode, and it is a phase" begin
    npar = 9
    base = _saem_prior_objective(npar)
    mk() = ctsem_laplace_objective(base; re_index=[1, 2, 5], sd_index=[6, 7, 9],
        cor_index=[8], sd_scale=[1.0, 1.0, 1.0], level_nre=[2, 1],
        group=vcat(1:6, _TWOLEVEL_GROUP), level_ngroups=[6, 3])
    start = [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25]
    plain = ctsem_laplace_optimize(mk(), start; maxiter=500, progress=false)
    mode = plain.minimizer
    H = ctsem_laplace_hessian(mk(), mode)
    se = sqrt.(diag(inv(Symmetric(-H))))
    run = ctsem_saem(mk(), start .+ 0.3; seed=2)
    # Settled on the trend rule, well inside the cap.
    @test run.settled
    @test run.iterations < 10000
    @test run.trend <= _S._saem_trend_null()
    # The average is within a third of a standard error of the exact mode in
    # every coordinate: SAEM's Monte Carlo error is small beside the
    # estimate's own.
    @info "SAEM against the exact mode" maximum(abs.(run.minimizer .- mode) ./ se)
    @test all(abs.(run.minimizer .- mode) .<= se ./ 3)
    # Same seed, same run; another seed, other draws.
    again = ctsem_saem(mk(), start .+ 0.3; maxiter=200, seed=2)
    again2 = ctsem_saem(mk(), start .+ 0.3; maxiter=200, seed=2)
    other = ctsem_saem(mk(), start .+ 0.3; maxiter=200, seed=3)
    @test again.minimizer == again2.minimizer
    @test again.minimizer != other.minimizer
    # As a phase of the optimiser: SAEM runs, then L-BFGS polishes to the same
    # mode the plain optimiser finds.
    phased = ctsem_laplace_optimize(mk(), start; maxiter=500, progress=false,
        saem=true, saem_maxiter=400, saem_seed=4)
    @test phased.saem_iterations > 0
    @test phased.saem_trace !== nothing
    @test isapprox(phased.minimizer, mode; atol=1e-3)
    @test isapprox(phased.maximum_loglik, plain.maximum_loglik; atol=1e-6)
    off = ctsem_laplace_optimize(mk(), start; maxiter=500, progress=false)
    @test off.saem_iterations == 0
end

@testset "SAEM: the split across workers does not change the run" begin
    laplace, values = _fresh_twolevel()
    previous = ctsem_max_chunks().max_chunks
    try
        runs = map((1, max(2, Threads.nthreads()))) do w
            ctsem_set_max_chunks!(w)
            l, v = _fresh_twolevel()
            ctsem_saem(l, v; maxiter=60, seed=9).minimizer
        end
        @test runs[1] == runs[2]
    finally
        ctsem_set_max_chunks!(previous)
    end
end

# Three subject-level effects at rank one: loadings at raw 6-8, no scales.
_saem_reduced() = (ctsem_laplace_objective(_saem_prior_objective(8);
    re_index=[1, 2, 5], sd_index=Int[], cor_index=Int[], sd_scale=[1.0, 1.0, 1.0],
    level_nre=[3], group=collect(1:6), level_ngroups=[6], level_rank=[1],
    load_index=[6, 7, 8]), [0.2, -0.1, 0.3, -0.2, 0.05, 0.6, 0.4, 0.3])

@testset "SAEM: fixed point under a reduced-rank level" begin
    laplace, values = _saem_reduced()
    plain = ctsem_laplace_optimize(laplace, values; maxiter=500, progress=false)
    run = ctsem_saem(_saem_reduced()[1], values; seed=6)
    @test run.settled
    # The covariance, not the loadings: a rank-one loading's sign is not
    # identified, and either sign is the same model.
    cov_of(x) = ctsem_laplace_popcov(laplace, x, 1)
    @test isapprox(cov_of(run.minimizer), cov_of(plain.minimizer); rtol=0.15, atol=0.02)
    # Every other parameter within a third of its standard error, one by one.
    se = sqrt.(diag(inv(Symmetric(-ctsem_laplace_hessian(laplace, plain.minimizer)))))
    @test all(abs.(run.minimizer[1:5] .- plain.minimizer[1:5]) .<= se[1:5] ./ 3)
end

@testset "SAEM: chains" begin
    @test _S._saem_default_chains(6) == 8
    @test _S._saem_default_chains(13) == 4
    @test _S._saem_default_chains(800) == 1
    laplace, values = _fresh_twolevel()          # three units: eight chains
    st = ctsem_saem_init(laplace, values; seed=4)
    @test st.chains == 8
    for _ in 1:20
        ctsem_saem_step!(st, laplace)
    end
    # Every chain started at its own draw; after twenty iterations no two
    # chains of a unit hold the same one.
    for U in eachindex(laplace.units.members)
        @test length(unique(st.u[U])) == st.chains
    end
    # One chain, asked for, is one chain.
    st1 = ctsem_saem_init(laplace, values; seed=4, chains=1)
    @test st1.chains == 1
end

# The centred step re-expresses the draws: theta's population means and
# scales move, and every member's shifted parameters stay exactly where they
# were, so no likelihood moves -- on one, two and three levels and at reduced
# rank. And with no prior on a level's means, the step brings that level's
# effects' average over its groups and chains toward zero, the Gaussian's
# maximum -- part of the way, by the data's share of each direction.
@testset "SAEM: the centred step leaves every member where it was" begin
    function shifted_all(st, laplace)
        Ls = _S._laplace_popchols(st.theta, laplace.spec)
        [_S._laplace_member_values(st.theta, laplace.spec, Ls, st.u[U][c],
            laplace.units.offsets[U][m])
         for U in eachindex(st.u) for c in 1:st.chains
         for m in eachindex(laplace.units.members[U])]
    end
    for (name, fresh) in (("one level", _fresh_linear), ("two levels", _fresh_twolevel),
                          ("three levels", _fresh_threelevel), ("rank one", _saem_reduced))
        laplace, values = fresh()
        st = ctsem_saem_init(laplace, values; seed=3)
        for _ in 1:15
            ctsem_saem_step!(st, laplace; centre=false)
        end
        pop = reduce(vcat, [_S._laplace_level_positions(laplace.spec, l)
                            for l in eachindex(laplace.spec.levels)])
        keep = setdiff(eachindex(st.theta), pop)
        prec = _S._saem_prior_precision(laplace, length(st.theta))
        before = shifted_all(st, laplace)
        theta0 = copy(st.theta)
        total_before = map(eachindex(laplace.spec.levels)) do l
            t = zeros(_S.nlatent(laplace.spec.levels[l]))
            for U in eachindex(st.u), blk in laplace.units.blocks[U], c in 1:st.chains
                blk.level == l || continue
                t .+= st.u[U][c][(blk.offset + 1):(blk.offset + blk.size)]
            end
            t
        end
        _S._saem_centre!(st, laplace, prec)
        after = shifted_all(st, laplace)
        @test st.theta != theta0
        @test all(isapprox(a[keep], b[keep]; atol=1e-12) for (a, b) in zip(after, before))
        for (l, level) in enumerate(laplace.spec.levels)
            all(iszero, prec[level.re_index]) || continue
            r = _S.nlatent(level)
            total = zeros(r)
            for U in eachindex(st.u), blk in laplace.units.blocks[U], c in 1:st.chains
                blk.level == l || continue
                total .+= st.u[U][c][(blk.offset + 1):(blk.offset + blk.size)]
            end
            @test norm(total) <= norm(total_before[l]) + 1e-10
        end
    end
end

# The reduced-rank loading step's factor: without a prior the closed form
# chol(S / G); with one, a stationary point of the same objective.
@testset "SAEM: the loading factor maximises its objective" begin
    rng = Random.Xoshiro(3)
    r, k, G = 2, 4, 30
    U = randn(rng, r, G) .* [1.3, 0.7]
    S = U * transpose(U)
    R = zeros(k, r)
    for q in 1:r, p in q:k
        R[p, q] = 0.5 + 0.1 * (p + q)
    end
    A0 = Matrix(cholesky(Symmetric(S ./ G)).L)
    @test _S._saem_loading_factor(A0, S, R, zeros(k, r), G) ≈ A0 atol = 1e-8
    D = zeros(k, r)
    for q in 1:r, p in q:k
        D[p, q] = 4.0
    end
    A = _S._saem_loading_factor(A0, S, R, D, G)
    idx = [(p, q) for q in 1:r for p in q:r]
    f(x) = begin
        B = zeros(eltype(x), r, r)
        for (j, (p, q)) in enumerate(idx)
            B[p, q] = x[j]
        end
        Binv = inv(LowerTriangular(B))
        -G * sum(a -> log(abs(a)), diag(B)) - tr(Binv * S * transpose(Binv)) / 2 -
            sum(D .* (R * B) .^ 2) / 2
    end
    x = [A[p, q] for (p, q) in idx]
    @test maximum(abs, ForwardDiff.gradient(f, x)) < 1e-6
    @test f(x) > f([A0[p, q] for (p, q) in idx])
end

# A loading pattern (`load_row`, `load_col`; `ctFit(popblocks =)`): a move
# `L -> L A` uses only the entries the pattern allows, so it maps the pattern
# to itself and leaves every deviation `L u` where it was.
@testset "SAEM: a loading pattern's moves keep to it" begin
    # Two blocks over four effects: effect 1 loads on both, 4 on the first,
    # 2 and 3 on the second.
    level = _S.CTSEMLaplaceLevel(1:4, Int[], Int[], ones(4), [1, 2], 2;
        rank=2, load_index=1:5, load_row=[1, 4, 1, 2, 3], load_col=[1, 1, 2, 2, 2])
    @test _S.isreducedrank(level)
    @test _S._laplace_loading_entries(level) == ([1, 4, 1, 2, 3], [1, 1, 2, 2, 2])
    allowed = _S._laplace_loading_allowed(level)
    @test allowed == [true false; false true]
    L = _S._laplace_poploading(collect(1.0:5.0), level)
    @test L == [1.0 3.0; 0.0 4.0; 0.0 5.0; 2.0 0.0]
    # The triangle allows its whole lower triangle.
    tri = _S.CTSEMLaplaceLevel(1:4, Int[], Int[], ones(4), [1, 2], 2;
        rank=2, load_index=1:7)
    @test _S._laplace_loading_allowed(tri) == [true false; true true]
    rng = Random.Xoshiro(4)
    M = randn(rng, 2, 40)
    S = M * transpose(M)
    A, whole = _S._saem_block_cholesky(S ./ 40, allowed)
    @test whole && A[1, 2] == 0 && A[2, 1] == 0
    @test A[1, 1] ≈ sqrt(S[1, 1] / 40) && A[2, 2] ≈ sqrt(S[2, 2] / 40)
    Atri, wtri = _S._saem_block_cholesky(S ./ 40, _S._laplace_loading_allowed(tri))
    @test wtri && Atri ≈ Matrix(cholesky(Symmetric(S ./ 40)).L)
    E = [0.3 0.0; 0.0 -0.2]
    LA = L * exp(E)
    @test LA[2, 1] == 0 && LA[3, 1] == 0 && LA[4, 2] == 0
    u = randn(rng, 2)
    @test LA * (exp(-E) * u) ≈ L * u
    # A coordinate is inert where every effect loading on it is.
    @test !_S._laplace_coordinate_inert(level, 2, [2, 3])
    @test _S._laplace_coordinate_inert(level, 2, [1, 2, 3])
    @test _S._laplace_coordinate_inert(level, 1, [1, 4])
    @test !_S._laplace_coordinate_inert(level, 1, [4])
end

# The SAEM sampler draws the same joint posterior NUTS draws (`ctsem_sample`,
# validated on its own): means within four Monte Carlo standard errors,
# standard deviations within ten per cent, at two draw counts, since a wrong
# sampler returns a plausible wrong answer whose error does not shrink. Two
# fixtures, so every move is in it: the rank-one level (the loadings'
# transformation move), and two nested full-rank levels (the centred and
# non-centred scale moves, the collapsed moves of a block above a leaf).
function _saem_against_nuts(mk, values; ndraws)
    mode = ctsem_laplace_optimize(mk(), values; maxiter=500, progress=false).minimizer
    nuts = ctsem_sample(mk(), mode; nchains=4, nwarmup=300, ndraws=1500, seed=3)
    map(ndraws) do n
        s = ctsem_saem_sample(mk(), mode; nchains=4, nwarmup=300, ndraws=n, seed=7)
        a = nuts.draws; b = s.draws
        sa = vec(std(a; dims=2)); sb = vec(std(b; dims=2))
        mcse = sqrt.(sa .^ 2 ./ nuts.ess .+ sb .^ 2 ./ s.ess)
        (z=maximum(abs.(vec(mean(b; dims=2)) .- vec(mean(a; dims=2))) ./ mcse),
         sd=extrema(sb ./ sa), rhat=maximum(s.rhat), sampler=s.sampler,
         ndraws=s.ndraws, ncp=s.ncp_accept)
    end
end

_saem_twolevel() = ctsem_laplace_objective(_saem_prior_objective(9);
    re_index=[1, 2, 5], sd_index=[6, 7, 9], cor_index=[8], sd_scale=[1.0, 1.0, 1.0],
    level_nre=[2, 1], group=vcat(1:6, _TWOLEVEL_GROUP), level_ngroups=[6, 3])

# The linear fixture with TI predictors: an estimated effect (raw 9) of the
# first on CINT, one of the random-effect parameters, and an effect of the
# second fixed at 0.5 on T0MEANS, the other. Ten subjects.
isdefined(@__MODULE__, :_SAEM_TI_OBJECTIVE) || global _SAEM_TI_OBJECTIVE = let
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
            (:MANIFESTMEANS, 1, 1) => (5, "param[5]")))
    ti = DataFrame(parameter=[2, 1], predictor=[1, 2], coefficient=[9, 0],
        value=[0.0, 0.5])
    starts, times, data = _laplace_test_data(10, 1)
    x = hcat(collect(range(-1.0, 1.0; length=10)), [sin(1.3 * s) for s in 1:10])
    _S.CTSEMObjective(ekf_from_data_frame(df, ti), starts, times, data,
        zeros(0, size(data, 2)), x; prior_index=collect(1:9), prior_scale=ones(9))
end
_saem_ti() = (ctsem_laplace_objective(_SAEM_TI_OBJECTIVE, [1, 2], [6, 7], [8],
    [1.0, 1.0]), [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, 0.3])

# Every chain's members' parameters as the filter sees them, TI effects added.
function _saem_ti_members(st, laplace; chains=st.chains)
    Ls = _S._laplace_popchols(st.theta, laplace.spec)
    out = Vector{Vector{Float64}}()
    for U in eachindex(st.u), c in 1:chains, (m, s) in enumerate(laplace.units.members[U])
        v = _S._laplace_member_values(st.theta, laplace.spec, Ls, st.u[U][c],
            laplace.units.offsets[U][m])
        so = laplace.objective.subject_objectives[s]
        push!(out, _S._materialize_subject_values!(Float64[], v, so.params, so.tipreds))
    end
    return out
end

# A TI coefficient on a random-effect parameter moves with the level's mean
# in the centred step and in the sampler's centred draw, and the effects take
# its change back out, so no member's parameters move; a fixed effect has no
# coefficient to move.
@testset "SAEM: TI coefficients move with the centred mean; no member moves" begin
    laplace, values = _saem_ti()
    @test _S._saem_level_ti(laplace, 1) == [(9, [(1, 2)])]
    st = ctsem_saem_init(laplace, values; seed=3)
    for _ in 1:15
        ctsem_saem_step!(st, laplace; centre=false)
    end
    prec = _S._saem_prior_precision(laplace, length(st.theta))
    model = 1:5
    before = _saem_ti_members(st, laplace)
    theta0 = copy(st.theta)
    _S._saem_centre!(st, laplace, prec)
    after = _saem_ti_members(st, laplace)
    @test st.theta[9] != theta0[9]
    @test all(isapprox(a[model], b[model]; atol=1e-12) for (a, b) in zip(after, before))
    before = _saem_ti_members(st, laplace; chains=1)
    theta0 = copy(st.theta)
    _S._saem_sample_mean!(st, laplace, 1, prec, Random.Xoshiro(5))
    after = _saem_ti_members(st, laplace; chains=1)
    @test st.theta[9] != theta0[9]
    @test all(isapprox(a[model], b[model]; atol=1e-12) for (a, b) in zip(after, before))
end

# Natural deviations (`d = S u`, what a sampled fit reports and generation takes)
# must rebuild exactly the subject parameters the sampler's own coordinates
# give -- nested, and on a reduced-rank level whose deviations are wider than
# its coordinates.
@testset "deviations rebuild the parameters the effects give" begin
    for (laplace, values) in ((_saem_twolevel(),
            [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25]), _saem_reduced())
        theta = collect(Float64, values)
        Ls = _S._laplace_popchols(theta, laplace.spec)
        units = laplace.units
        effects = randn(Random.Xoshiro(4), sum(units.dims))
        starts, ndev = _S._laplace_deviation_offsets(laplace)
        dev = _S._laplace_deviations!(zeros(ndev), laplace, Ls, effects, starts)
        ustarts = cumsum(units.dims) .- units.dims
        fromu = _S._laplace_shifted_values(laplace, theta, Ls, U ->
            effects[(ustarts[U] + 1):(ustarts[U] + units.dims[U])])
        fromd = _S._laplace_shifted_values_from_deviations(laplace, theta, dev)
        @test fromd ≈ fromu
        layout = _S.ctsem_laplace_deviation_layout(laplace)
        @test length(layout.position) == ndev
        @test all(>(0), layout.within) && all(>(0), layout.unit)
    end
end

# The samplers standardise each level's effects in a Cholesky ordered
# strongest-determined first (`_with_root_orders`). Any order gives the same
# covariance, Laplace value and gradient -- det(I + S'HS) = det(I + H Sigma) --
# and the natural order is restored after a run.
@testset "ordered roots: same covariance, value and gradient; cleared after a run" begin
    laplace = _saem_twolevel()
    values = [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25]
    base = ctsem_laplace_evaluate(laplace, values; gradient=true)
    cov0 = [L * transpose(L) for L in _S._laplace_popchols(values, laplace.spec)]
    orders = _S._laplace_order_roots!(laplace, values)
    shares = _S._laplace_determined_shares(laplace, values)
    @test all(0 .<= shares[1] .<= 1)
    @test orders[1] == sortperm(shares[1]; rev=true)
    @test isempty(orders[2])                     # one effect: nothing to order
    copyto!(laplace.spec.levels[1].order, [2, 1])  # a permuted root, whatever the shares
    Ls = _S._laplace_popchols(values, laplace.spec)
    @test !istril(Ls[1])
    @test Ls[1] * transpose(Ls[1]) ≈ cov0[1]
    ev = ctsem_laplace_evaluate(laplace, values; gradient=true)
    @test ev.value ≈ base.value rtol = 1e-8
    @test ev.gradient ≈ base.gradient rtol = 1e-6
    _S._laplace_clear_root_orders!(laplace)
    @test all(isempty(level.order) for level in laplace.spec.levels)
    s = ctsem_saem_sample(laplace, values; nchains=2, nwarmup=30, ndraws=30, seed=3)
    @test all(isfinite, s.draws)
    @test all(isempty(level.order) for level in laplace.spec.levels)
end

@testset "SAEM sampler: the joint posterior NUTS draws ($name)" for (name, mk, values) in (
        ("rank one", () -> _saem_reduced()[1], _saem_reduced()[2]),
        ("two full-rank levels", _saem_twolevel,
         [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25]),
        ("TI coefficient on a random-effect parameter", () -> _saem_ti()[1],
         _saem_ti()[2]))
    r = _saem_against_nuts(mk, values; ndraws=(1000, 4000))
    @info "SAEM sampler against NUTS ($name)" r
    @test all(x -> x.sampler == "saem", r)
    # No target: exactly the draws asked for.
    @test [x.ndraws for x in r] == [1000, 4000]
    @test r[2].z < 4
    @test 0.9 < r[2].sd[1] && r[2].sd[2] < 1.1
    @test r[2].rhat < 1.05
    # No discrepancy that persists as the draws grow.
    @test r[2].z < max(4, r[1].z)
    name == "two full-rank levels" && @test all(isfinite, r[2].ncp)
end

# The same stopping rule as NUTS (`_sample_until_target`): a target the first
# batch cannot meet extends the run, within the budget, and the same seed
# reproduces it at one worker. Above one the pool divides each evaluation's
# sum over units as workers come free, so its last bits differ from run to
# run and the chains, which amplify them, part (CLAUDE.md, "A julia test that
# fails on one machine..."); this asserted reproducibility at any width and
# failed whenever the suite ran with threads.
@testset "SAEM sampler: sampled to a target, and reproducible" begin
    laplace, values = _saem_reduced()
    mode = ctsem_laplace_optimize(laplace, values; maxiter=500, progress=false).minimizer
    held = ContinuousTimeSEM.ctsem_max_chunks().max_chunks
    ContinuousTimeSEM.ctsem_set_max_chunks!(1)
    s, again = try
        (ctsem_saem_sample(_saem_reduced()[1], mode; nchains=2, nwarmup=100,
                ndraws=40, min_ess=150, max_draws=600, seed=11),
            ctsem_saem_sample(_saem_reduced()[1], mode; nchains=2, nwarmup=100,
                ndraws=40, min_ess=150, max_draws=600, seed=11))
    finally
        ContinuousTimeSEM.ctsem_set_max_chunks!(held)
    end
    @test 40 < s.ndraws <= 600
    @test size(s.draws, 2) == 2 * s.ndraws
    @test s.min_ess >= 150 || s.ndraws == 600
    @test again.draws == s.draws
    # The effects' summaries are over the same draws, as natural deviations:
    # wider than the rank-one coordinates, one per parameter the level moves.
    sampler = ctsem_sampler(_saem_reduced()[1], length(mode))
    ndev = last(_S._laplace_deviation_offsets(_saem_reduced()[1]))
    @test ndev == 6 * 3 && ndev > sampler.ndim - sampler.npar
    @test length(s.effect_mean) == ndev
    @test all(isfinite, s.effect_sd)
    saved = ctsem_saem_sample(_saem_reduced()[1], mode; nchains=2, nwarmup=50,
        ndraws=30, seed=11, save_effects=true)
    @test size(saved.draws, 1) == sampler.npar + ndev
end

# Every sampler takes its chains' starts from a caller -- one column a chain --
# and reports when each check of the stopping rule happened; the SAEM kernel
# also continues a SAEM state's chains rather than starting fresh.
@testset "samplers take starts, the SAEM kernel continues a SAEM state, each reports its checks" begin
    laplace, values = _saem_reduced()
    run = ctsem_saem(laplace, values; seed=2)
    npar = length(values)
    starts = repeat(run.minimizer, 1, 2) .+ 0.01 .* randn(Random.Xoshiro(1), npar, 2)

    s = ctsem_saem_sample(_saem_reduced()[1], run.minimizer; nchains=2, nwarmup=60,
        ndraws=40, seed=3, state=run.state, min_ess=50, max_draws=200)
    @test s.sampler == "saem"
    tr = s.target_trace
    @test length(tr.draws) >= 1
    @test issorted(tr.secs) && all(>=(0), tr.secs)
    @test tr.draws[end] == s.ndraws
    @test all(isfinite, s.draws)

    withstarts = ctsem_saem_sample(_saem_reduced()[1], run.minimizer; nchains=2,
        nwarmup=20, ndraws=20, seed=3, starts=starts)
    @test size(withstarts.draws) == (npar, 40)
    @test isempty(withstarts.target_trace.draws)          # no target, no checks
    @test_throws DimensionMismatch ctsem_saem_sample(_saem_reduced()[1], run.minimizer;
        nchains=2, nwarmup=5, ndraws=5, starts=zeros(3, 2))

    nuts = ctsem_sample(_saem_reduced()[1], run.minimizer; nchains=2, nwarmup=20,
        ndraws=20, seed=3, starts=starts, min_ess=30, max_draws=100)
    @test length(nuts.target_trace.draws) >= 1
    @test nuts.target_trace.draws[end] == nuts.ndraws
    marginal = ctsem_sample_marginal(_saem_reduced()[1], run.minimizer; nchains=2,
        nwarmup=20, ndraws=20, seed=3, starts=starts)
    @test size(marginal.draws) == (npar, 40)
    @test all(isfinite, marginal.draws)
end

# `saem = true`: both joint samplers place themselves from SAEM's own run --
# its estimate, its chains' effects -- rather than from the point handed in,
# which need only be a start.
@testset "both joint samplers place from SAEM's state when asked" begin
    laplace, values = _saem_reduced()
    s = ctsem_saem_sample(laplace, values; nchains=2, nwarmup=40, ndraws=30, seed=4,
        saem=true)
    p = s.placement
    @test p !== nothing
    @test p.saem_iterations > 0
    @test length(p.theta) == length(values) && all(isfinite, p.theta)
    @test all(isfinite, s.draws)
    n = ctsem_sample(_saem_reduced()[1], values; nchains=2, nwarmup=40, ndraws=30,
        seed=4, saem=true)
    @test n.placement !== nothing && n.placement.saem_iterations > 0
    @test all(isfinite, n.draws)
    # Without it, nothing is placed and nothing is reported.
    @test ctsem_saem_sample(_saem_reduced()[1], values; nchains=2, nwarmup=5,
        ndraws=5, seed=4).placement === nothing
    # `init_scale` scales the chains' jitter about SAEM's estimate, as it
    # scales NUTS's own starts when nothing is placed: the same run, so the
    # same estimate and the same draw, twice as far out at 2, and every chain
    # at the estimate at 0. It was accepted and ignored here before.
    p1 = _S._saem_placement(_saem_reduced()[1], values; nchains=2, seed=4)
    p2 = _S._saem_placement(_saem_reduced()[1], values; nchains=2, seed=4,
        init_scale=2.0)
    p0 = _S._saem_placement(_saem_reduced()[1], values; nchains=2, seed=4,
        init_scale=0.0)
    @test p2.theta == p1.theta
    @test any(p1.pars .!= p1.theta)
    @test p2.pars .- p2.theta ≈ 2 .* (p1.pars .- p1.theta)
    @test all(p0.pars .== p0.theta)
    run = ctsem_saem(_saem_reduced()[1], values; seed=2, maxiter=200)
    @test_throws ArgumentError ctsem_saem_sample(_saem_reduced()[1], values; nchains=2,
        nwarmup=5, ndraws=5, saem=true, state=run.state)
end
