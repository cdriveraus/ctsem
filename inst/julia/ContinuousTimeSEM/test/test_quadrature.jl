using ForwardDiff, LinearAlgebra

# `_fresh_linear`, `_fresh_nonlinear`, `_fresh_twolevel` and `_fresh_threelevel`
# come from `laplace_fixtures.jl`, which `test_laplace.jl` also includes. Shared
# rather than rebuilt because `ekf_from_columns` parses transform strings through
# `eval`, so a rebuilt model is uniquely typed and recompiles the whole filter.
isdefined(@__MODULE__, :_LAPLACE_LINEAR_OBJECTIVE) ||
    include(joinpath(@__DIR__, "laplace_fixtures.jl"))

# Adaptive Gauss-Hermite quadrature over the same integral `laplace.jl`
# approximates.
#
# Three things need checking, and they are independent of one another.
#
#  1. The *rule* is a Gauss-Hermite rule. Checked against integrals whose value
#     is known in closed form, with no ctsem in sight.
#
#  2. The rule and the Laplace term are integrating the *same thing*. Checked at
#     `nodes = 1`, where the adaptive rule degenerates to a single point at the
#     mode and must reproduce the Laplace value exactly -- constants, log
#     determinant, prior and all. That is a sharp test: the two arrive at the
#     number by completely different routes, and every constant that either one
#     drops has to be dropped by the other.
#
#  3. Where Laplace is exact, more nodes must change nothing. The linear model's
#     integrand is exactly Gaussian in `z`, so a nine-node rule integrates the
#     same Gaussian and has to agree to quadrature precision.

@testset "Gauss-Hermite rule" begin
    x, w = ContinuousTimeSEM._gauss_hermite(7)
    @test length(x) == 7 && length(w) == 7
    @test sum(w) ≈ sqrt(pi) rtol = 1e-12
    # int x^2 exp(-x^2) dx = sqrt(pi)/2, and the rule is exact to degree 2m-1.
    @test sum(w .* x .^ 2) ≈ sqrt(pi) / 2 rtol = 1e-12
    @test sum(w .* x .^ 4) ≈ 3 * sqrt(pi) / 4 rtol = 1e-12
    @test sum(w .* x .^ 3) ≈ 0 atol = 1e-12
    # One node is the degenerate rule the adaptive scheme uses as Laplace.
    x1, w1 = ContinuousTimeSEM._gauss_hermite(1)
    @test x1 == [0.0] && w1 ≈ [sqrt(pi)]

    points, logweights = ContinuousTimeSEM._gh_grid(2, 5)
    @test length(points) == 25 && length(logweights) == 25
    # The exp(x^2) factor is folded into the log weight, so the plain weights
    # come back by removing it.
    @test sum(exp(lw - sum(p .^ 2)) for (p, lw) in zip(points, logweights)) ≈ pi rtol = 1e-12
end

@testset "quadrature reproduces Laplace at one node" begin
    for (name, fresh) in (("linear", _fresh_linear), ("nonlinear", _fresh_nonlinear))
        laplace, values = fresh()
        reference = ctsem_laplace_evaluate(laplace, values; gradient=false).value
        single = ctsem_laplace_quadrature(laplace, values; nodes=1).value
        @test single ≈ reference rtol = 1e-12
        @test isfinite(single)
        name == "linear" || continue
        # Where Laplace is exact the rule cannot move: the integrand is exactly
        # Gaussian in z, so every node count integrates the same Gaussian.
        for nodes in (3, 5, 9)
            @test ctsem_laplace_quadrature(laplace, values; nodes=nodes).value ≈
                reference rtol = 1e-9
        end
    end
end

@testset "quadrature disagrees with Laplace only where it should" begin
    laplace, values = _fresh_nonlinear()
    reference = ctsem_laplace_evaluate(laplace, values; gradient=false).value
    five = ctsem_laplace_quadrature(laplace, values; nodes=5).value
    nine = ctsem_laplace_quadrature(laplace, values; nodes=9).value
    # A genuinely non-Gaussian integrand: more nodes must move the answer, and
    # the movement must be settling rather than wandering.
    @test five != reference
    @test abs(nine - five) < abs(five - reference)
    @test isapprox(nine, ctsem_laplace_quadrature(laplace, values; nodes=11).value;
        atol = 1e-6)
end

@testset "quadrature contributions sum to the value" begin
    laplace, values = _fresh_nonlinear()
    result = ctsem_laplace_quadrature(laplace, values; nodes=5, contributions=true)
    prior = ContinuousTimeSEM._ctsem_log_prior(laplace.objective, collect(Float64, values))
    @test sum(result.subject_loglik) + prior ≈ result.value rtol = 1e-12
    @test length(result.subject_loglik) == length(laplace.objective.subject_objectives)
end

@testset "the block tree the recursion walks" begin
    for (label, fresh, nblocks, nroots) in (("two levels", _fresh_twolevel, 3, 1),
                                            ("three levels", _fresh_threelevel, 7, 1))
        laplace, _ = fresh()
        blocks = laplace.units.blocks[1]
        tree = ContinuousTimeSEM._quadrature_children(blocks)
        @test (label, length(blocks)) == (label, nblocks)
        @test (label, length(tree.roots)) == (label, nroots)
        # Every block is either a root or exactly one block's child, and the
        # tree spans all of them -- otherwise the recursion would silently drop
        # a level's likelihood contribution.
        reached = Set(tree.roots)
        queue = copy(tree.roots)
        while !isempty(queue)
            b = pop!(queue)
            for c in tree.children[b]
                @test !(c in reached)
                push!(reached, c); push!(queue, c)
            end
        end
        @test (label, length(reached)) == (label, length(blocks))
        # A child sits one level below its parent, and its members are a subset.
        for b in eachindex(blocks), c in tree.children[b]
            @test blocks[c].level < blocks[b].level
            @test issubset(Set(blocks[c].members), Set(blocks[b].members))
        end
        # A leaf's members partition its ancestors' members.
        leaves = [b for b in eachindex(blocks) if isempty(tree.children[b])]
        covered = sort(vcat([blocks[b].members for b in leaves]...))
        @test covered == collect(1:length(laplace.units.members[1]))
    end
end

@testset "the recursion reproduces Laplace at one node, at every depth" begin
    for (label, fresh) in (("two levels", _fresh_twolevel),
                           ("three levels", _fresh_threelevel))
        laplace, values = fresh()
        reference = ctsem_laplace_evaluate(laplace, values; gradient=false).value
        single = ctsem_laplace_quadrature(laplace, values; nodes=1).value
        @test (label, isapprox(single, reference; rtol=1e-10)) == (label, true)
        # `nodes = 1` forces `sum over blocks of logdet(scale) = -logdet(M)/2`,
        # and that sum telescopes only if each outer block is scaled by its
        # curvature *after* its descendants are eliminated. The marginal
        # covariance passes this at two levels and fails at three, which is why
        # the check runs at both depths rather than one.
    end
    # Two levels of *identity*-transformed effects make the integrand exactly
    # Gaussian in the whole unit latent vector, so every node count integrates
    # the same Gaussian -- through the recursion, the per-node leaf mode solves
    # and all. The three-level fixture puts its outermost effect on DRIFT, whose
    # transform is not the identity, so it is genuinely approximate and belongs
    # to the test below instead.
    laplace, values = _fresh_twolevel()
    reference = ctsem_laplace_evaluate(laplace, values; gradient=false).value
    for nodes in (3, 5)
        @test ctsem_laplace_quadrature(laplace, values; nodes=nodes).value ≈
            reference rtol = 1e-7
    end
end

@testset "three levels: a nonlinear outer effect moves the answer and settles" begin
    laplace, values = _fresh_threelevel()
    reference = ctsem_laplace_evaluate(laplace, values; gradient=false).value
    five = ctsem_laplace_quadrature(laplace, values; nodes=5).value
    nine = ctsem_laplace_quadrature(laplace, values; nodes=9).value
    # The outermost effect sits on DRIFT, so Laplace is approximate here and the
    # rule must disagree with it -- and then stop disagreeing with itself.
    @test five != reference
    @test abs(nine - five) < abs(five - reference)
    @test isapprox(nine, ctsem_laplace_quadrature(laplace, values; nodes=11).value;
        atol = 1e-5)
end

@testset "per-unit contributions at more than one level" begin
    laplace, values = _fresh_twolevel()
    result = ctsem_laplace_quadrature(laplace, values; nodes=3, contributions=true)
    prior = ContinuousTimeSEM._ctsem_log_prior(laplace.objective, collect(Float64, values))
    # One entry per *unit*, not per subject: a study's integral does not
    # decompose over its members.
    @test length(result.subject_loglik) == length(laplace.units.members)
    @test length(result.subject_loglik) < length(laplace.objective.subject_objectives)
    @test sum(result.subject_loglik) + prior ≈ result.value rtol = 1e-12
end

@testset "the correction is available at more than one level" begin
    laplace, values = _fresh_twolevel()
    result = ctsem_laplace_optimize(laplace, values; maxiter=100, g_tol=1e-7,
        tune_chunks=false)
    H = ctsem_laplace_hessian(laplace, result.minimizer)
    correction = ctsem_laplace_correction(laplace, result.minimizer, H; nodes=3)
    @test all(isfinite, correction.delta)

    # Not `-H * delta == gap` any more, and the reason is the sentence below
    # about `inv(-H)` being "big enough to carry it": on this fixture the gap
    # gradient is entirely floating-point noise (Laplace is exact here) and the
    # information matrix is near-singular, so an exact solve amplifies the noise
    # by 1/lambda_min. Two consecutive corrections on the same fit returned
    # largest steps of 3.3e+08 and 3.5e-10 standard errors. `_correction_step`
    # truncates at sqrt(eps) relative to the largest eigenvalue, so the residual
    # now lives in the directions it declined to correct along.
    #
    # Three properties replace the one, and together they are stronger.

    # 1. Where the curvature *is* identified, the solve is still exact --
    #    truncation must be a no-op on a well-conditioned problem, or it has
    #    quietly broken every model that did not need it.
    n = length(result.minimizer)
    wellposed = ctsem_laplace_correction(laplace, result.minimizer,
        Matrix(-1.0I, n, n); nodes=3)
    @test wellposed.dropped_directions == 0
    @test -(Symmetric(Matrix(-1.0I, n, n)) * wellposed.delta) ≈
        wellposed.gap_gradient rtol = 1e-7

    # 2. A singular information matrix yields no correction rather than an
    #    arbitrarily large one.
    singular = ctsem_laplace_correction(laplace, result.minimizer,
        zeros(n, n); nodes=3)
    @test all(iszero, singular.delta)
    @test singular.dropped_directions == n

    # 3. The correction on the real Hessian is reproducible. This is what
    #    actually failed before: the two calls disagreed by eighteen orders of
    #    magnitude while describing the same fit.
    again = ctsem_laplace_correction(laplace, result.minimizer, H; nodes=3)
    @test correction.delta ≈ again.delta atol = 1e-9
    # Laplace is exact for this model, so there is nothing to correct. The step
    # is measured against the standard errors rather than against zero: the gap
    # it differences is ~1e-9 here, a 1e-3 finite-difference step turns that
    # into ~5e-7 of gradient, and six subjects make `inv(-H)` big enough to
    # carry it to 1e-3 in raw units -- which is still four orders of magnitude
    # inside the interval, and it is the ratio that a reader acts on.
    se = sqrt.(abs.(diag(inv(-Symmetric(H)))))
    @test maximum(abs.(correction.delta ./ se)) < 0.05
    @test abs(correction.gap) < 1e-6
end

@testset "the correction is the Newton step it claims to be" begin
    laplace, values = _fresh_nonlinear()
    result = ctsem_laplace_optimize(laplace, values; maxiter=200, g_tol=1e-8,
        tune_chunks=false)
    est = result.minimizer
    H = ctsem_laplace_hessian(laplace, est)
    correction = ctsem_laplace_correction(laplace, est, H; nodes=5)
    @test all(isfinite, correction.delta)
    @test correction.corrected ≈ est .+ correction.delta
    # `delta` solves -H delta = grad(Q - T). At the Laplace estimate the Laplace
    # gradient is zero, so that gap gradient is the quadrature objective's own
    # gradient, and the step is an ordinary Newton step on it.
    @test -(Symmetric(H) * correction.delta) ≈ correction.gap_gradient rtol = 1e-8
    # The gap gradient must be a derivative of the gap the same call reports.
    @test correction.gap ≈ correction.quadrature - correction.laplace rtol = 1e-12

    # Refinement moves in the direction the correction predicts, and improves
    # the objective it is maximising.
    refined = ctsem_laplace_refine(laplace, est; nodes=5, maxiter=30)
    @test refined.maximum_loglik >=
        ctsem_laplace_quadrature(laplace, est; nodes=5).value - 1e-8
    @test refined.shift ≈ refined.minimizer .- est
end

@testset "quadrature is unaffected by how the modes were left" begin
    # The rule reuses the retained inner modes as its centre, so it inherits
    # whatever `_laplace_solve_unit_mode!` guarantees. That guarantee is the
    # point: a detour through a distant parameter vector must not change the
    # value here.
    laplace, values = _fresh_nonlinear()
    first = ctsem_laplace_quadrature(laplace, values; nodes=5).value
    detour = collect(Float64, values); detour[2] += 3.0; detour[5] += 2.5
    ctsem_laplace_quadrature(laplace, detour; nodes=5)
    @test ctsem_laplace_quadrature(laplace, values; nodes=5).value ≈ first rtol = 1e-10
end
