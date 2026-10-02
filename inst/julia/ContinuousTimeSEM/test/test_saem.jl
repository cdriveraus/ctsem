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
function _saem_estep_draws(laplace, values, U, n; seed=1)
    st = ctsem_saem_init(laplace, values; seed=seed)
    Ls = _S._laplace_popchols(st.theta, laplace.spec)
    d = length(st.u[U][1])
    draws = zeros(d, n)
    for k in 1:n
        st.iteration += 1
        _S._laplace_parallel(laplace, [U]) do U
            (k % 25 == 0) && _S._saem_refresh!(st, laplace, U, Ls)
            _S._saem_sweep!(st, laplace, U, 1, Ls, 1, 1 / (1 + k)^0.6)
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

@testset "SAEM: E-step draws the exact conditional ($name)" for (name, fresh) in (
        ("one level", _fresh_linear), ("two levels", _fresh_twolevel),
        ("three levels", _fresh_threelevel))
    laplace, values = fresh()
    U = 1
    mode, Sigma = _saem_exact_conditional(laplace, values, U)
    sd = sqrt.(diag(Sigma))
    errs = map((4000, 16000)) do n
        draws, st = _saem_estep_draws(laplace, values, U, n + 500; seed=7)
        x = draws[:, 501:end]
        m = vec(mean(x; dims=2))
        C = cov(x; dims=2)
        (mean=maximum(abs.(m .- mode) ./ sd),
         var=maximum(abs.(diag(C) ./ diag(Sigma) .- 1)),
         cor=maximum(abs.(cov2cor(C) .- cov2cor(Sigma))),
         acc=_S._saem_acceptance(st))
    end
    @info "SAEM E-step moments ($name)" errs
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
    # Settled on the drift rule, well inside the cap, with Kesten's steps
    # smaller at the end than at the start and never growing.
    @test run.settled
    @test run.iterations < 10000
    @test run.trace.gamma[1] == 1
    @test all(diff(run.trace.gamma) .<= 1e-12)
    @test run.trace.gamma[end] < 1
    @test run.drift < 0.01
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
    @test isapprox(run.minimizer[1:5], plain.minimizer[1:5]; atol=0.1)
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
    # Every chain started at the mode; after twenty iterations no two chains
    # of a unit hold the same draw.
    for U in eachindex(laplace.units.members)
        @test length(unique(st.u[U])) == st.chains
    end
    # One chain, asked for, is one chain.
    st1 = ctsem_saem_init(laplace, values; seed=4, chains=1)
    @test st1.chains == 1
end
