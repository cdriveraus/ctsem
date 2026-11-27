using ForwardDiff, LinearAlgebra, Random, Statistics

# The variational objective (`variational.jl`): each unit's effects integrated
# by a lower bound, its approximation fitted at every evaluation.
#
# What has to hold, each checked on its own:
#
#  1. On a Gaussian integrand the Gaussian family's bound is the integral,
#     which Laplace gives exactly there -- on units of one block, of subjects
#     in studies and of three levels. The tree family is exact there with
#     parameters linear in the members, not quadratic in the coordinates.
#  2. The shaped family's log Jacobian is log |det du/dz| of its own map, its
#     hand-written second derivatives are the generic map's, and its entropy
#     on its own rule is the nodes' on a rule fine enough to settle it.
#  3. The gradient is the derivative of the value (the envelope theorem), on
#     non-Gaussian integrands: one effect a subject, and subjects in studies.
#  4. On a non-Gaussian integrand the bound lies below the integral, which
#     brute-force quadrature gives.
#  5. The family chosen by need, and an optimisation, its value after the run
#     the run's.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

const _VF = ContinuousTimeSEM

# Subjects nested in two studies on the nonlinear objective (as test_placed.jl's
# `_fresh_placed_nested`): a T0MEANS effect a subject, a MANIFESTMEANS effect a
# study.
_vf_nested_nonlinear() = (ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE;
    re_index=[1, 4], sd_index=[5, 6], cor_index=Int[], sd_scale=[1.0, 1.0],
    level_nre=[1, 1], group=vcat(1:5, [1, 1, 2, 2, 2]), level_ngroups=[5, 2]),
    [0.1, 0.4, -0.2, 0.05, -0.25, -0.1])

# Three levels with every effect a mean's (the fixtures' `_fresh_threelevel`
# puts the region's on the drift, through a nonlinear transform): subjects'
# T0MEANS and CINT, studies' and regions' MANIFESTMEANS.
_vf_threelevel_gauss() = (ctsem_laplace_objective(_LAPLACE_LINEAR_OBJECTIVE;
    re_index=[1, 2, 5, 5], sd_index=[6, 7, 9, 10], cor_index=[8],
    sd_scale=[1.0, 1.0, 1.0, 1.0], level_nre=[2, 1, 1],
    group=vcat(1:6, _THREELEVEL_STUDY, _THREELEVEL_REGION), level_ngroups=[6, 3, 2]),
    [0.2, -0.1, 0.3, -0.2, 0.05, -0.3, -0.15, 0.4, -0.25, 0.1])

@testset "variational objective" begin
    @testset "Gaussian integrand: the bound is the integral" begin
        for fresh in (_fresh_linear, _fresh_twolevel, _vf_threelevel_gauss)
            laplace, values = fresh()
            exact = ctsem_laplace_evaluate(laplace, values; gradient=false)
            V = ctsem_variational(laplace; family="gauss", tol=1e-12)
            r = _VF.ctsem_variational_evaluate(V, values)
            @test r.converged
            @test r.unit_loglik ≈ exact.unit_loglik atol = 1e-8
        end
        # A skewed or bending family on a Gaussian integrand stays at zero
        # shape up to its rule's error.
        laplace, values = _fresh_linear()
        exact = ctsem_laplace_evaluate(laplace, values; gradient=false)
        for family in ("skew", "bend")
            V = ctsem_variational(laplace; family=family, tol=1e-12)
            r = _VF.ctsem_variational_evaluate(V, values)
            @test r.converged
            @test maximum(abs.(r.unit_loglik .- exact.unit_loglik)) < 0.003
        end
    end

    @testset "the tree family's parameters are linear in the members" begin
        for fresh in (_fresh_twolevel, _fresh_threelevel)
            laplace, _ = fresh()
            V = ctsem_variational(laplace; family="gauss")
            for (U, lay) in enumerate(V.layouts)
                blocks = laplace.units.blocks[U]
                want = sum(b.size * (b.size + 3) ÷ 2 +
                           b.size * sum((blocks[a].size for a in b.ancestors); init=0)
                           for b in blocks; init=0)
                @test lay.p == want
                @test length(lay.pblocks) == length(blocks)
                # Each member's term reads its own chain's parameters alone.
                for mem in lay.members
                    mine = [b for b in eachindex(blocks) if mem.member in blocks[b].members]
                    @test all(q -> lay.pblock[q] in mine, mem.plist)
                    @test length(mem.chain) == sum(blocks[b].size for b in mine)
                end
            end
        end
    end

    @testset "shaped family: log Jacobian, second derivatives, entropy" begin
        rng = Random.Xoshiro(3)
        for (d, skew, bend) in ((1, true, false), (2, true, false), (2, true, true),
            (3, true, true), (3, false, true))
            f = _VF.CTSEMVFamily(d, skew; bend=bend)
            lam = 0.3 .* randn(rng, f.p)
            Z = randn(rng, d, 5)
            nd = _VF._vf_nodes(f, lam, Z)
            for k in 1:5
                z = Z[:, k]
                Jz = ForwardDiff.jacobian(zz -> _VF._vf_map_generic(f, lam, zz)[1], z)
                @test nd.lj[k] ≈ logabsdet(Jz)[1] atol = 1e-10
                u, l = _VF._vf_map_generic(f, lam, z)
                @test nd.U[:, k] ≈ u atol = 1e-12
                @test nd.J[k] ≈ ForwardDiff.jacobian(_VF._VFMapAt(f, z), lam) atol = 1e-10
            end
            Zf, wf = _VF._vf_rule_product(d, d == 1 ? 200 : d == 2 ? 100 : 30)
            @test _VF._vf_entropy(f, lam) ≈
                  dot(wf, _VF._vf_nodes(f, lam, Zf; jacobian=false).lj) atol = 2e-5
            if skew && !bend
                r = randn(rng, d)
                H = zeros(f.p, f.p)
                _VF._vf_second!(H, 1.0, f, lam, _VF._vf_L(f, lam), Z[:, 1], r)
                @test (H .+ transpose(H)) ./ 2 ≈
                      ForwardDiff.hessian(_VF._VFUAt(f, Z[:, 1], r), lam) atol = 1e-8
            end
        end
    end

    @testset "the gradient is the value's" begin
        for (fresh, family) in ((_fresh_nonlinear, "auto"), (_fresh_nonlinear, "gauss"),
            (_vf_nested_nonlinear, "auto"))
            laplace, values = fresh()
            V = ctsem_variational(laplace; family=family, tol=1e-14)
            r = _VF.ctsem_variational_evaluate(V, values)
            @test r.converged
            h = 1e-4
            for j in eachindex(values)
                up = copy(values)
                up[j] += h
                dn = copy(values)
                dn[j] -= h
                fd = (_VF.ctsem_variational_evaluate(V, up; gradient=false).value -
                      _VF.ctsem_variational_evaluate(V, dn; gradient=false).value) / (2h)
                @test r.gradient[j] ≈ fd rtol = 1e-4 atol = 1e-6
            end
        end
    end

    @testset "below the integral" begin
        laplace, values = _fresh_nonlinear()
        V = ctsem_variational(laplace)
        r = _VF.ctsem_variational_evaluate(V, values)
        z, w = _VF._vf_gh(80)
        for U in eachindex(laplace.units.members)
            a = [log(w[k]) + _unit_loglik(laplace, U, values, [z[k]]) for k in eachindex(z)]
            pk = maximum(a)
            exact = pk + log(sum(exp.(a .- pk)))
            @test r.unit_loglik[U] <= exact + 1e-4
            # A skewed approximation of one effect leaves little of the integral.
            @test exact - r.unit_loglik[U] < 0.01
        end
    end

    @testset "the family chosen by need" begin
        # One block of two coordinates: started Gaussian, tried with the
        # bending family where a run stops. On a Gaussian integrand bending
        # gains nothing, so no unit switches and each stays a candidate; on a
        # non-Gaussian one the bound cannot fall.
        laplace, values = _fresh_linear()
        V = ctsem_variational(laplace)
        @test all(V.pending) && all(l -> l.shaped === nothing, V.layouts)
        _VF.ctsem_variational_evaluate(V, values)
        c = ctsem_variational_choose!(V, values)
        @test c.units == 0 && all(V.pending)
        # The nonlinear fixture with effects on T0MEANS and MANIFESTMEANS:
        # raw 1-4 the model, 5-6 the scales, 7 their correlation.
        laplace2 = ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE, [1, 4], [5, 6], [7],
            [1.0, 1.0])
        values2 = [0.1, 0.4, -0.2, 0.05, -0.25, -0.4, 0.2]
        V2 = ctsem_variational(laplace2; need=0.0)
        r2 = _VF.ctsem_variational_evaluate(V2, values2)
        c2 = ctsem_variational_choose!(V2, values2)
        r3 = _VF.ctsem_variational_evaluate(V2, values2)
        @test c2.units == count(l -> l.shaped !== nothing && l.shaped.bend, V2.layouts)
        @test r3.value >= r2.value - 1e-8
        @test r3.value ≈ r2.value + c2.change atol = 1e-6
        ctsem_variational_reset!(V2)
        @test all(V2.pending) && all(l -> l.shaped === nothing, V2.layouts)
        # A nested unit's family is the tree's and is not tried.
        V3 = ctsem_variational(_vf_nested_nonlinear()[1])
        @test !any(V3.pending)
    end

    @testset "the Hessian is the solved bound's" begin
        # The implicit Hessian against differences of gradients solved to the
        # last bit, where no unit's approximation stays fixed at the points; at
        # the fixtures' values, since the formula holds wherever each unit's
        # approximation is at its optimum. Then the solves started from the
        # accepted approximation moved by the Hessian's derivative reach the
        # same optimum as solves from the accepted approximation alone.
        laplace2 = ctsem_laplace_objective(_LAPLACE_NONLINEAR_OBJECTIVE, [1, 4], [5, 6], [7],
            [1.0, 1.0])
        values2 = [0.1, 0.4, -0.2, 0.05, -0.25, -0.4, 0.2]
        for (laplace, values, family) in ((_fresh_nonlinear()..., "auto"),
            (_vf_nested_nonlinear()..., "auto"), (laplace2, values2, "gauss"),
            (laplace2, values2, "bend"))
            V = ctsem_variational(laplace; family=family, tol=1e-14)
            _VF.ctsem_variational_evaluate(V, values)
            Hd = _VF.ctsem_hessian(V, values; method=:differenced)
            Hi = _VF.ctsem_hessian(V, values)
            @test Hi ≈ Hd rtol = 1e-4 atol = 1e-5
            @test V.dlam !== nothing
            moved = values .+ 0.05 .* (-1) .^ (1:length(values))
            rp = _VF.ctsem_variational_evaluate(V, moved; gradient=false)
            W = ctsem_variational(laplace; family=family, tol=1e-14)
            _VF.ctsem_variational_evaluate(W, values)
            ra = _VF.ctsem_variational_evaluate(W, moved; gradient=false)
            @test rp.unit_loglik ≈ ra.unit_loglik atol = 1e-8
        end
    end

    @testset "the probe's screen samples the units' bounds" begin
        # A coreset's units alone, at its weights, leaving the trial; every
        # unit's weighted gradient as the scores; and the screen's trial
        # value the full one where it takes every unit.
        laplace, values = _fresh_nonlinear()
        n = length(values)
        nunits = length(laplace.units.members)
        V = ctsem_variational(laplace; tol=1e-12)
        full = _VF.ctsem_variational_evaluate(V, values)
        prior = full.value - sum(full.unit_loglik)
        trial = deepcopy(V.lam_trial)
        cs = _VF._laplace_coreset(V, [2, 4], [1.0, 3.0])
        sub = _VF._laplace_coreset_try(V, values, false, cs)
        @test sub.unit_loglik[[2, 4]] ≈ full.unit_loglik[[2, 4]] atol = 1e-10
        @test all(isnan, sub.unit_loglik[setdiff(1:nunits, [2, 4])])
        @test sub.value ≈ full.unit_loglik[2] + 3 * full.unit_loglik[4] + prior atol = 1e-10
        @test V.lam_trial == trial
        S = zeros(nunits, n)
        r = _VF._laplace_coreset_try(V, values, true, nothing; scores=S)
        gp = zeros(length(_VF._laplace_extend(laplace, values)))
        _VF._ctsem_log_prior_gradient!(gp, laplace.objective, _VF._laplace_extend(laplace, values))
        @test vec(sum(S; dims=1)) .+ gp[1:n] ≈ r.gradient atol = 1e-10
        @test V.lam_trial == trial
        s = _VF.CTSEMProbeScreen(V, values, nothing, Float64[], Float64[], NaN, NaN,
            falses(0, 0), zeros(0, 0))
        xt = values .+ 0.05 .* (-1.0) .^ (1:n)
        @test _VF._ctsem_probe_trial(s, V, xt, full.value, Inf) ≈
              _VF.ctsem_variational_evaluate(V, xt; gradient=false).value atol = 1e-8
    end

    @testset "a stage taken back" begin
        # A run that R does not keep leaves its approximations accepted; the
        # snapshot taken before it puts back the kept stage's, so the
        # estimate's value needs no solve.
        laplace, values = _fresh_nonlinear()
        V = ctsem_variational(laplace)
        res = ctsem_optimize(V, values; maxiter=200, newton=false)
        x = collect(Float64, res.minimizer)
        before = _VF.ctsem_variational_evaluate(V, x)
        snap = _VF.ctsem_variational_snapshot(V)
        ctsem_optimize(V, values .+ 0.3; maxiter=5, newton=false)
        @test V.theta_acc != x
        _VF.ctsem_variational_restore!(V, snap)
        @test V.theta_acc == snap.state.theta_acc
        i0 = _VF.ctsem_variational_info(V)
        after = _VF.ctsem_variational_evaluate(V, x)
        i1 = _VF.ctsem_variational_info(V)
        @test after.value ≈ before.value atol = 1e-10
        @test after.gradient ≈ before.gradient atol = 1e-8
        @test i1.held_steps == i0.held_steps
        @test i1.exact_iterations == i0.exact_iterations
    end

    @testset "optimised, and read after the run" begin
        for fresh in (_fresh_nonlinear, _vf_nested_nonlinear)
            laplace, values = fresh()
            V = ctsem_variational(laplace)
            res = ctsem_optimize(V, values; maxiter=200, newton=false)
            x = collect(Float64, res.minimizer)
            after = _VF.ctsem_variational_evaluate(V, x; gradient=true)
            @test after.value ≈ res.maximum_loglik atol = 1e-6
            @test norm(after.gradient) < 1e-3
            @test V.accepted && !V.running
        end
    end

    @testset "units' terms and rows' scores read the bound" begin
        # On a Gaussian integrand the bound is Laplace's term, so leave-one-out's
        # units' terms and the score tests' row scores are Laplace's; and
        # neither moves the approximations the objective holds.
        laplace, values = _fresh_linear()
        V = ctsem_variational(laplace; family="gauss", tol=1e-12)
        _VF.ctsem_variational_evaluate(V, values)
        held = deepcopy(V.lam)
        draws = hcat(values, values .+ 0.05)
        tv = ctsem_unit_terms(V, draws)
        tl = ctsem_unit_terms(laplace, draws)
        @test all(tv.converged)
        @test tv.unit_loglik ≈ tl.unit_loglik atol = 1e-7
        @test tv.value ≈ tl.value atol = 1e-6
        one = ctsem_unit_terms(V, draws[:, 1:1]; subjects=2)
        @test one.unit_loglik[2, 1] ≈ tl.unit_loglik[2, 1] atol = 1e-7
        @test all(isnan, one.unit_loglik[[1; 3:end], 1]) && isnan(one.value[1])
        n = length(values)
        D = Matrix{Float64}(I, n, n)
        @test ctsem_row_scores(V, values, D) ≈ ctsem_row_scores(laplace, values, D) atol = 1e-6
        @test V.lam == held
        # Off it a subject's rows sum to its unit's score on the bound.
        laplace, values = _fresh_nonlinear()
        V = ctsem_variational(laplace; tol=1e-12)
        r = _VF.ctsem_variational_evaluate(V, values; units=true)
        D = Matrix{Float64}(I, length(values), length(values))
        R = ctsem_row_scores(V, values, D)
        @test all(isfinite, R)
        rows = [size(s.data, 2) for s in laplace.objective.full_subjects]
        stop = cumsum(rows)
        for (U, members) in enumerate(laplace.units.members)
            i = members[1]
            span = (stop[i] - rows[i] + 1):stop[i]
            w = _VF._laplace_unit_weight(laplace, U)
            @test vec(sum(R[span, :]; dims=1)) .* w ≈ r.unit_gradient[:, U] atol = 1e-6
        end
        # The Hessian's columns are the whole Hessian's, and leave the
        # approximations' derivatives where the whole one put them.
        Hfull = ctsem_hessian(V, values)
        dlam = deepcopy(V.dlam)
        Hc = ctsem_hessian_columns(V, values, [1, 3])
        @test Hc ≈ Hfull[:, [1, 3]] rtol = 1e-4 atol = 1e-5
        @test V.dlam == dlam
    end

    @testset "ctLaplaceCheck's gap is against the bound" begin
        # The quadrature's integral less the bound (here within the
        # quadrature's own error of zero).
        laplace, values = _fresh_nonlinear()
        V = ctsem_variational(laplace; tol=1e-12)
        _VF.ctsem_variational_evaluate(V, values)
        H = ctsem_hessian(V, values)
        c = ctsem_laplace_correction(V, values, H; nodes=9)
        @test c.gap ≈ c.quadrature - c.laplace
        @test c.laplace ≈ _VF.ctsem_variational_evaluate(V, values; gradient=false,
            keep=false).value atol = 1e-10
        @test abs(c.gap) < 0.01
    end
end
