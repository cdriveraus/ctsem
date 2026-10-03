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
        ones(npar), 1.0, Int[], Float64[], Float64[])
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

# The centred move's conditional of a full-rank level's scales and correlations
# given the deviations, whose gradient is written by hand: against central
# differences, on the two-effect level (scales and the correlation) and the
# one-effect level above it.
@testset "SAEM sampler: the centred move's gradient" begin
    laplace = _saem_twolevel()
    values = [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25]
    st = _S.ctsem_saem_init(laplace, values; seed=5, chains=1)
    prec = _S._saem_prior_precision(laplace, length(values))
    for l in eachindex(laplace.spec.levels)
        positions, _, _, density! = _S._saem_centred_density(st, laplace, l, prec)
        x = st.theta[positions] .+ 0.1
        g = zeros(length(x))
        @test isfinite(density!(g, x))
        h = 1e-6
        fd = map(eachindex(x)) do t
            e = zeros(length(x)); e[t] = h
            (density!(zeros(length(x)), x .+ e) - density!(zeros(length(x)), x .- e)) / (2h)
        end
        @test g ≈ fd rtol = 1e-5 atol = 1e-6
    end
end

@testset "SAEM sampler: the joint posterior NUTS draws ($name)" for (name, mk, values) in (
        ("rank one", () -> _saem_reduced()[1], _saem_reduced()[2]),
        ("two full-rank levels", _saem_twolevel,
         [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25]))
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
# reproduces it whatever the thread count did.
@testset "SAEM sampler: sampled to a target, and reproducible" begin
    laplace, values = _saem_reduced()
    mode = ctsem_laplace_optimize(laplace, values; maxiter=500, progress=false).minimizer
    s = ctsem_saem_sample(_saem_reduced()[1], mode; nchains=2, nwarmup=100,
        ndraws=40, min_ess=150, max_draws=600, seed=11)
    @test 40 < s.ndraws <= 600
    @test size(s.draws, 2) == 2 * s.ndraws
    @test s.min_ess >= 150 || s.ndraws == 600
    again = ctsem_saem_sample(_saem_reduced()[1], mode; nchains=2, nwarmup=100,
        ndraws=40, min_ess=150, max_draws=600, seed=11)
    @test again.draws == s.draws
    # The effects' summaries are over the same draws, in the joint layout.
    sampler = ctsem_sampler(_saem_reduced()[1], length(mode))
    @test length(s.effect_mean) == sampler.ndim - sampler.npar
    @test all(isfinite, s.effect_sd)
    saved = ctsem_saem_sample(_saem_reduced()[1], mode; nchains=2, nwarmup=50,
        ndraws=30, seed=11, save_effects=true)
    @test size(saved.draws, 1) == sampler.ndim
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
