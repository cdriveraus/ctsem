# Discrete-time models.
#
# Discrete time is not a second filter here, it is a different discretization:
# a discrete model's DRIFT, CINT and DIFFUSION are already the one-step
# quantities, so the matrix exponential, the Lyapunov solve and the intercept
# solve that turn continuous parameters into per-interval ones all collapse.
# Everything downstream -- the measurement update, the smoother, generation --
# is untouched.
#
# So the tests are: the forward pass really is the plain one-step recursion, the
# reverse pass agrees with forward-mode differentiation of it (which is the only
# independent check of a hand-written adjoint), and the recorded times play no
# part, which is the property that distinguishes discrete from continuous.

function _discrete_setup(; times=[0.0, 1.0, 2.0, 3.0, 4.0], continuous=false)
    cells = [
        (:T0MEANS, 1, 1, 1, missing, "param[1]", missing, missing),
        (:LAMBDA, 1, 1, missing, 1.0, missing, missing, missing),
        (:DRIFT, 1, 1, 2, missing, "1 / (1 + exp(-param[2]))", missing, missing),
        (:DIFFUSION, 1, 1, 3, missing, "log1p_exp(param[3])", missing, missing),
        (:MANIFESTVAR, 1, 1, 4, missing, "log1p_exp(param[4])", missing, missing),
        (:MANIFESTMEANS, 1, 1, 5, missing, "param[5]", missing, missing),
        (:CINT, 1, 1, 6, missing, "param[6]", missing, missing),
        (:T0VAR, 1, 1, 7, missing, "log1p_exp(param[7])", missing, missing),
        (:JAx, 1, 1, missing, missing, missing, "DRIFT[1,1]", missing),
        (:Jy, 1, 1, missing, 1.0, missing, missing, missing),
        (:PARS, 1, 1, missing, 0.0, missing, missing, missing),
    ]
    df = DataFrame(
        matrix = [c[1] for c in cells],
        row = [c[2] for c in cells],
        col = [c[3] for c in cells],
        parnumber = Union{Missing,Int}[c[4] for c in cells],
        value = Union{Missing,Float64}[c[5] for c in cells],
        transform = Union{Missing,String}[c[6] for c in cells],
        predicttransform = Union{Missing,String}[c[7] for c in cells],
        updatetransform = Union{Missing,String}[c[8] for c in cells],
    )
    sp = ContinuousTimeSEM.ekf_from_columns(
        df.matrix, df.row, df.col,
        [ismissing(v) ? 0 : Int(v) for v in df.parnumber],
        [ismissing(v) ? NaN : Float64(v) for v in df.value],
        [ismissing(v) ? "" : String(v) for v in df.transform],
        [ismissing(v) ? "" : String(v) for v in df.predicttransform],
        [ismissing(v) ? "" : String(v) for v in df.updatetransform],
        fill("", nrow(df));
        continuous_time=continuous)
    data = reshape([0.4, -0.2, 0.5, NaN, 0.1], 1, 5)
    objective = ctsem_objective(sp, [1], times, data)
    return (sp=sp, objective=objective, times=times,
        values=[0.3, 0.4, -0.1, -0.6, 0.2, 0.05, -0.3])
end

_discrete_plain = _discrete_setup()
_discrete_stretched = _discrete_setup(times=[0.0, 2.5, 7.0, 9.0, 20.0])
_continuous_same = _discrete_setup(continuous=true)

@testset "the discrete forward pass is the plain one-step recursion" begin
    setup = _discrete_plain
    k = ctsem_kalman(setup.objective, setup.values)
    layout = ctsem_parameter_layout(setup.objective)
    flat = ctsem_parameter_matrices(setup.objective, setup.values)
    block(name) = begin
        j = findfirst(==(name), layout.matrix)
        reshape(flat[(layout.offset[j]+1):(layout.offset[j]+layout.nrow[j]*layout.ncol[j]), 1],
            layout.nrow[j], layout.ncol[j])
    end
    A = block("DRIFT")
    Q = block("DIFFUSIONcov")
    CINT = block("CINT")

    # x_prior = A x_upd + CINT, exactly -- no exponential, no discrete intercept
    # solve.
    for r in 2:5
        @test k.eta[1, r, :] ≈ A * k.eta[2, r-1, :] .+ CINT[:, 1] rtol = 1e-12
        # P_prior = A (P_upd + eps I) A' + Q, with the same 1e-10 ridge the
        # continuous path bakes into every propagation.
        Pupd = k.etacov[2, r-1, :, :] .+ 1e-10 .* [1.0;;]
        @test k.etacov[1, r, :, :] ≈ A * Pupd * A' .+ Q rtol = 1e-8
    end

    # The transition handed to the smoother is JAx itself.
    JAx = block("JAx")
    for r in 2:5
        @test k.transition[r, :, :] ≈ JAx rtol = 1e-12
    end

    # The asymptotic forms are the discrete ones: (I - A) x = c and
    # X = A X A' + Q.
    asymC = block("asymCINT")
    asymD = block("asymDIFFUSIONcov")
    @test (I - A) * asymC ≈ CINT rtol = 1e-10
    @test asymD ≈ A * asymD * A' .+ Q rtol = 1e-10
end

@testset "recorded times do not matter in discrete time" begin
    # The defining difference from a continuous model: each row advances exactly
    # one step whatever the interval says. Stretching the times must change
    # nothing at all.
    plain = ctsem_evaluate(_discrete_plain.objective, _discrete_plain.values;
        gradient=true)
    stretched = ctsem_evaluate(_discrete_stretched.objective,
        _discrete_stretched.values; gradient=true)
    @test plain.value ≈ stretched.value rtol = 1e-14
    @test plain.gradient ≈ stretched.gradient rtol = 1e-12

    # ... whereas the same model read as continuous does depend on them, so the
    # comparison above is not vacuous.
    continuous = ctsem_evaluate(_continuous_same.objective, _continuous_same.values;
        gradient=false)
    @test !isapprox(continuous.value, plain.value; rtol=1e-6)
end

# F5: `_reverse_predict_discrete!` (`adjoint_ekf.jl:1005-1060`) is a second
# implementation of the prediction reverse, and `_discrete_plain` /
# `_discrete_stretched` above are 1-latent, where every matrix is 1x1 and its
# four transposition-sensitive lines (`:1015,1017,1023,1024`) are all the
# identity. These reuse the 2-latent scenario builders from
# `test_adjoint_gradient_validation.jl` (non-diagonal DRIFT, and for the
# second, free off-diagonal DIFFUSION/MANIFESTVAR/T0VAR) with
# `continuous=false`, so a transposition or `Ps A P̃` swap in the discrete
# branch would show up here where it cannot at one latent state.
_discrete_cross_2d = (
    objective=ctsem_objective(
        _adjoint_cross_effect_2d_parameters(continuous=false), [1],
        [0.0, 0.5, 1.2], reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2], 2, :)),
    values=[0.2, 0.1, -0.3],
)
_discrete_free_covariance = (
    objective=ctsem_objective(
        _adjoint_free_covariance_2d_parameters(continuous=false), [1],
        [0.0, 0.5, 1.2, 2.0],
        reshape([0.1, -0.2, 0.15, 0.05, -0.1, 0.2, 0.3, -0.05], 2, :)),
    values=[0.3, 0.2, -0.1, 0.4, 0.15, -0.2, 0.25, 0.1, -0.15, 0.35],
)

@testset "the discrete adjoint agrees with forward-mode differentiation" begin
    # The only independent check of a hand-written reverse pass. ForwardDiff
    # differentiates the primal itself, so an error in the discrete branch of
    # the adjoint cannot hide in both.
    for setup in (_discrete_plain, _discrete_stretched, _discrete_cross_2d,
        _discrete_free_covariance)
        adjoint = ctsem_adjoint_gradient(setup.objective, setup.values)
        forward = ctsem_evaluate(setup.objective, setup.values; gradient=true,
            gradient_method="forward")
        @test adjoint.value ≈ forward.value rtol = 1e-12
        @test adjoint.gradient ≈ forward.gradient rtol = 1e-8
    end

    # And per-subject scores still sum to it, so the discrete branch did not
    # break the accumulation the score path depends on.
    multi = ctsem_objective(_discrete_plain.sp, [1, 4], [0.0, 1.0, 2.0, 0.0, 1.0],
        reshape([0.4, -0.2, 0.5, 0.1, -0.3], 1, 5))
    summed = ctsem_adjoint_gradient(multi, _discrete_plain.values)
    per_subject = ctsem_subject_gradients(multi, _discrete_plain.values)
    @test vec(sum(per_subject.scores, dims=1)) ≈ summed.gradient rtol = 1e-10
end

@testset "prediction and generation work in discrete time" begin
    setup = _discrete_plain
    k = ctsem_kalman(setup.objective, setup.values)
    # Smoothing still runs backwards over the subject, and the last row is the
    # filtered estimate.
    @test k.eta[3, 5, :] ≈ k.eta[2, 5, :]
    @test tr(k.etacov[3, 1, :, :]) < tr(k.etacov[2, 1, :, :])
    # The missing row updates nothing.
    @test k.eta[2, 4, :] ≈ k.eta[1, 4, :]

    g = ctsem_generate(setup.objective, setup.values,
        reshape([0.7, -1.2, 0.3, 0.9, 1.5], 1, 5))
    @test isnan(g.Y[1, 4])
    regenerated = ctsem_objective(setup.sp, [1], setup.times, g.Y)
    @test ctsem_evaluate(regenerated, setup.values; gradient=false).value ≈
        sum(g.subject_loglik) rtol = 1e-12
end
