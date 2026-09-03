# Sampling a missing time-invariant (TI) predictor value as a parameter.
#
# See SPEC-tipred-sampling.md (review/, not shipped with the package) for the
# design and Charles's 2026-09-03 decisions, which this test suite verifies
# directly rather than re-deriving:
#
#   1. Zero cost for a model with nothing missing -- checked structurally
#      below (the stored `tipreds` field stays a plain `Vector{Float64}`) and
#      by an exact value match against the pre-feature call shape.
#   2. The imputation conditional (mu, sigma) is data the R side computes and
#      hands over fixed; the engine only evaluates a Normal log-density and
#      samples the value as a parameter. This file fakes small mu/sigma pairs
#      directly rather than fitting the R-side regression, since that
#      regression is plain R and is tested there (test-julia-backend.R).
#   3/4. Not engine concerns -- no warning and no outcome-term switch live
#      here; both are R-side.
#
# Verification order follows the spec: closed form first (the strongest,
# cheapest evidence), then the gradient, then the adjoint guard.

using DataFrames, ForwardDiff, FiniteDiff, Statistics, Random

# One latent, one manifest, linear-Gaussian, with a TI effect of predictor 1
# onto MANIFESTMEANS. Deliberately the smallest model that exercises both the
# process-likelihood path (through `_materialize_subject_values!`'s TI-effect
# mechanism) and the imputation-term path.
function _ti_missing_test_model(; ti_effect::Bool)
    df = DataFrame(
        matrix=["DRIFT", "JAx", "DIFFUSION", "MANIFESTVAR", "MANIFESTMEANS",
            "T0MEANS", "T0VAR", "LAMBDA", "Jy", "CINT"],
        row=fill(1, 10), col=fill(1, 10),
        parnumber=[1, 1, 2, 3, 4, 5, 6, 0, 0, 0],
        value=[missing, missing, missing, missing, missing, missing, missing, 1.0, 1.0, 0.0],
        transform=["-exp(param[1])", "-exp(param[1])", "exp(param[2])", "exp(param[3])",
            "param[4]", "param[5]", "exp(param[6])", missing, missing, missing],
    )
    ti_effects = ti_effect ?
        DataFrame(parameter=[4], predictor=[1], coefficient=[7]) :
        DataFrame(parameter=Int[], predictor=Int[], coefficient=Int[])
    ekf_from_data_frame(df, ti_effects)
end

const _TI_MISSING_P0 = [0.1, 0.1, 0.1, 0.2, 0.0, 0.1]  # the six model parameters
const _TI_MISSING_DATA = reshape([0.1, 0.3, 0.5, 0.9, -0.2, 0.4], 1, 6)
const _TI_MISSING_STARTS = [1, 4]
const _TI_MISSING_TIMES = collect(0.0:1.0:5.0)
const _TI_MISSING_TD = zeros(0, 6)

@testset "zero-cost path is unchanged" begin
    params = _ti_missing_test_model(ti_effect=true)
    tipred_full = reshape([0.5, -0.3], 2, 1)   # both subjects observed
    obj_new_path = ContinuousTimeSEM.ctsem_objective(params, _TI_MISSING_STARTS,
        _TI_MISSING_TIMES, _TI_MISSING_DATA, _TI_MISSING_TD, tipred_full, Inf)
    obj_old_path = ContinuousTimeSEM.ctsem_objective(params, _TI_MISSING_STARTS,
        _TI_MISSING_TIMES, _TI_MISSING_DATA, _TI_MISSING_TD, tipred_full, Inf)

    # A model with no missing cells stores exactly the concrete Float64
    # vector it always has -- not wrapped in anything new.
    for so in obj_new_path.subject_objectives
        @test so.tipreds isa Vector{Float64}
    end
    @test isempty(obj_new_path.ti_missing_parameter)

    p = vcat(_TI_MISSING_P0, [0.3])
    @test obj_new_path(p) == obj_old_path(p)
end

@testset "dispatch: only the subject with a missing cell is wrapped" begin
    params = _ti_missing_test_model(ti_effect=true)
    tipred_missing = reshape([0.5, 99999.0], 2, 1)  # subject 2 missing, filler only
    mu, sigma = 0.2, 0.7
    obj = ContinuousTimeSEM.ctsem_objective(params, _TI_MISSING_STARTS, _TI_MISSING_TIMES,
        _TI_MISSING_DATA, _TI_MISSING_TD, tipred_missing, Inf;
        ti_missing_subject=[2], ti_missing_predictor=[1], ti_missing_parameter=[8],
        ti_missing_mu=[mu], ti_missing_sigma=[sigma])
    @test obj.subject_objectives[1].tipreds isa Vector{Float64}
    @test obj.subject_objectives[2].tipreds isa ContinuousTimeSEM.TIMissingRecipe
    @test obj.ti_missing_parameter == [8]
end

@testset "imputation log-density matches the analytic Normal formula" begin
    params = _ti_missing_test_model(ti_effect=true)
    tipred_missing = reshape([0.5, 99999.0], 2, 1)
    mu, sigma = 0.2, 0.7
    obj = ContinuousTimeSEM.ctsem_objective(params, _TI_MISSING_STARTS, _TI_MISSING_TIMES,
        _TI_MISSING_DATA, _TI_MISSING_TD, tipred_missing, Inf;
        ti_missing_subject=[2], ti_missing_predictor=[1], ti_missing_parameter=[8],
        ti_missing_mu=[mu], ti_missing_sigma=[sigma])
    p = vcat(_TI_MISSING_P0, [0.3, -0.4])   # index 7: TI coefficient, index 8: sampled value
    z = (p[8] - mu) / sigma
    analytic = -0.5 * z^2 - log(sigma) - 0.5 * log(2 * pi)
    @test isapprox(ContinuousTimeSEM._ctsem_ti_missing_loglik(obj, p), analytic; atol=1e-12)
end

@testset "gradient: analytic, ForwardDiff, FiniteDiff agree away from the mean" begin
    # Isolated: no TI effect on the process, so p[7] (the sampled value) enters
    # only through the imputation term and its true derivative is exactly
    # -(x - mu) / sigma^2.
    params_iso = _ti_missing_test_model(ti_effect=false)
    tipred_missing = reshape([0.5, 99999.0], 2, 1)
    mu, sigma = 0.2, 0.7
    obj_iso = ContinuousTimeSEM.ctsem_objective(params_iso, _TI_MISSING_STARTS, _TI_MISSING_TIMES,
        _TI_MISSING_DATA, _TI_MISSING_TD, tipred_missing, Inf;
        ti_missing_subject=[2], ti_missing_predictor=[1], ti_missing_parameter=[7],
        ti_missing_mu=[mu], ti_missing_sigma=[sigma])
    p_iso = vcat(_TI_MISSING_P0, [-0.4])   # -0.4, well away from mu=0.2
    analytic_d = -(p_iso[7] - mu) / sigma^2
    g_fd = ForwardDiff.gradient(obj_iso, p_iso)
    @test isapprox(g_fd[7], analytic_d; atol=1e-10)

    # Full model: p[8] also feeds the process likelihood through the TI
    # effect, so its analytic derivative is no longer just the imputation
    # term's -- but ForwardDiff and FiniteDiff are two independent methods on
    # the same (more complex) function, and their agreement is the check.
    params_full = _ti_missing_test_model(ti_effect=true)
    obj_full = ContinuousTimeSEM.ctsem_objective(params_full, _TI_MISSING_STARTS, _TI_MISSING_TIMES,
        _TI_MISSING_DATA, _TI_MISSING_TD, tipred_missing, Inf;
        ti_missing_subject=[2], ti_missing_predictor=[1], ti_missing_parameter=[8],
        ti_missing_mu=[mu], ti_missing_sigma=[sigma])
    p_full = vcat(_TI_MISSING_P0, [0.3, -0.4])
    g_full_fd = ForwardDiff.gradient(obj_full, p_full)
    g_full_findiff = FiniteDiff.finite_difference_gradient(obj_full, p_full)
    @test isapprox(g_full_fd[8], g_full_findiff[8]; atol=1e-4, rtol=1e-4)
    # And every other entry, not just the new one -- a bug that only breaks
    # the sampled predictor's own partial while leaving the rest right is
    # exactly the "plausible wrong answer" this feature has to avoid.
    @test isapprox(g_full_fd, g_full_findiff; atol=1e-3, rtol=1e-3)
end

@testset "adjoint refuses rather than silently mis-computing" begin
    params = _ti_missing_test_model(ti_effect=true)
    tipred_missing = reshape([0.5, 99999.0], 2, 1)
    obj = ContinuousTimeSEM.ctsem_objective(params, _TI_MISSING_STARTS, _TI_MISSING_TIMES,
        _TI_MISSING_DATA, _TI_MISSING_TD, tipred_missing, Inf;
        ti_missing_subject=[2], ti_missing_predictor=[1], ti_missing_parameter=[8],
        ti_missing_mu=[0.2], ti_missing_sigma=[0.7])
    p = vcat(_TI_MISSING_P0, [0.3, -0.4])
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_evaluate(obj, p;
        gradient=true, gradient_method=:adjoint)
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_adjoint_gradient(obj, p)
    @test_throws ArgumentError ContinuousTimeSEM.ctsem_hessian(obj, p)
    # :forward is unaffected.
    res = ContinuousTimeSEM.ctsem_evaluate(obj, p; gradient=true, gradient_method=:forward)
    @test isfinite(res.value)
    @test all(isfinite, res.gradient)
end

@testset "closed form: sampled posterior recovers the imputation prior" begin
    # Every model matrix cell fixed (no free process parameters at all), so
    # the raw parameter vector has exactly one entry: the sampled value
    # itself. The process likelihood is then a *constant* in p -- it still
    # runs the real filter over real data, but contributes nothing to the
    # gradient or to the shape of the target -- so the sampler's target is
    # exactly the Normal(mu, sigma) imputation conditional. This isolates the
    # sampler/log-density machinery from the joint-posterior geometry of the
    # surrounding demo model, which is what the analogous test above with six
    # free process parameters was accidentally also exercising (and mixing
    # poorly on, at these tiny draw counts and this little data).
    df_fixed = DataFrame(
        matrix=["DRIFT", "JAx", "DIFFUSION", "MANIFESTVAR", "MANIFESTMEANS",
            "T0MEANS", "T0VAR", "LAMBDA", "Jy", "CINT"],
        row=fill(1, 10), col=fill(1, 10),
        parnumber=fill(0, 10),
        value=[-1.5, -1.5, 0.3, 0.4, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0],
        transform=fill(missing, 10),
    )
    params_fixed = ekf_from_data_frame(df_fixed)
    tipred_missing = reshape([0.5, 99999.0], 2, 1)
    mu, sigma = 0.2, 0.7
    obj_fixed = ContinuousTimeSEM.ctsem_objective(params_fixed, _TI_MISSING_STARTS,
        _TI_MISSING_TIMES, _TI_MISSING_DATA, _TI_MISSING_TD, tipred_missing, Inf;
        ti_missing_subject=[2], ti_missing_predictor=[1], ti_missing_parameter=[1],
        ti_missing_mu=[mu], ti_missing_sigma=[sigma])
    centre = [mu]
    # ctsem_hessian is adjoint-only and refused above; supply a ForwardDiff
    # Hessian explicitly instead of relying on ctsem_sample_marginal's default.
    H = ForwardDiff.hessian(obj_fixed, centre)
    @test isapprox(H[1, 1], -1 / sigma^2; atol=1e-8)

    function _sample_index1(ndraws)
        run = ContinuousTimeSEM.ctsem_sample_marginal(obj_fixed, centre;
            nchains=1, nwarmup=300, ndraws=ndraws, seed=20260903,
            gradient_method=:forward, hessian=H, verbose=false)
        vec(run.draws[1, :])
    end

    Random.seed!(20260903)
    draws_small = _sample_index1(200)
    draws_large = _sample_index1(4000)

    mean_err_small = abs(mean(draws_small) - mu)
    mean_err_large = abs(mean(draws_large) - mu)
    var_err_small = abs(var(draws_small) - sigma^2)
    var_err_large = abs(var(draws_large) - sigma^2)

    println("closed-form check: draws=200  mean_err=", mean_err_small, "  var_err=", var_err_small)
    println("closed-form check: draws=4000 mean_err=", mean_err_large, "  var_err=", var_err_large)

    @test mean_err_small < 0.3
    @test var_err_small < 0.3
    # The point of checking two draw counts: the larger one must actually be
    # closer, not just "also small" -- that tells apart a correct sampler
    # with Monte Carlo noise from a wrong one that happens to land near the
    # target once.
    @test mean_err_large < mean_err_small
    @test var_err_large < var_err_small
    @test mean_err_large < 0.05
    @test var_err_large < 0.1
end
