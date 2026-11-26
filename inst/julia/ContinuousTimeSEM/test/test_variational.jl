using ForwardDiff, LinearAlgebra, Random, Statistics

# The variational objective (`variational.jl`): each unit's effects integrated
# by a lower bound, its approximation fitted at every evaluation.
#
# What has to hold, each checked on its own:
#
#  1. On a Gaussian integrand every family's bound is the integral, which
#     Laplace gives exactly there.
#  2. Each family's log Jacobian is log |det du/dz| of its own map, and the
#     skewed family's hand-written second derivatives are the generic map's.
#  3. The gradient is the derivative of the value (the envelope theorem: the
#     approximation maximises the bound, so its derivative in the parameters
#     is the partial at the nodes held), on a non-Gaussian integrand.
#  4. On a non-Gaussian integrand the bound lies below the integral, which
#     importance sampling on the fitted approximation gives.
#  5. The optimiser climbs it, and a value after the run is the run's.

isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

const _VF = ContinuousTimeSEM

@testset "variational objective" begin
    @testset "Gaussian integrand: the bound is the integral" begin
        laplace, values = _fresh_linear()
        exact = ctsem_laplace_evaluate(laplace, values; gradient=false)
        V = ctsem_variational(laplace; family="gauss", tol=1e-12)
        r = _VF.ctsem_variational_evaluate(V, values)
        @test r.converged
        @test r.unit_loglik ≈ exact.unit_loglik atol = 1e-8
        # A skewed or bending family on the cut rule leaves zero shape by the
        # rule's error, about a thousandth of a nat a unit here.
        for family in ("skew", "bend")
            V = ctsem_variational(laplace; family=family, tol=1e-12)
            r = _VF.ctsem_variational_evaluate(V, values)
            @test r.converged
            @test maximum(abs.(r.unit_loglik .- exact.unit_loglik)) < 0.003
        end
    end

    @testset "families: log Jacobian and second derivatives" begin
        rng = Random.Xoshiro(3)
        for (d, skew, bend) in ((1, true, false), (2, true, false), (2, true, true),
            (3, true, true), (3, false, true), (4, false, false))
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
            # The entropy on its fine rule against the nodes' log Jacobians
            # on a product rule fine enough to settle it.
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
        laplace, values = _fresh_nonlinear()
        for family in ("skew", "gauss")
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
        v = ctsem_variational_value(V, values; draws=8192)
        @test all(isfinite, v.unit_loglik)
        @test all(r.unit_loglik .<= v.unit_loglik .+ 4 .* v.unit_se .+ 1e-8)
        # A skewed approximation of one effect leaves little of the integral.
        @test maximum(v.unit_loglik .- r.unit_loglik) < 0.05
    end

    @testset "optimised, and read after the run" begin
        laplace, values = _fresh_nonlinear()
        V = ctsem_variational(laplace)
        res = ctsem_optimize(V, values; maxiter=200, newton=false)
        x = collect(Float64, res.minimizer)
        after = _VF.ctsem_variational_evaluate(V, x; gradient=true)
        @test after.value ≈ res.maximum_loglik atol = 1e-6
        @test norm(after.gradient) < 1e-3
        @test V.accepted && !V.running
    end
end
