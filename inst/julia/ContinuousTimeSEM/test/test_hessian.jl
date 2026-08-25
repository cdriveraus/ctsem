using DataFrames, ForwardDiff, LinearAlgebra

# The Hessian, by forward-mode differentiation of the reverse-mode gradient.
#
# `ForwardDiff.hessian` of the primal is the referee: it shares the primal with
# `ctsem_hessian` but nothing else -- no tape, no hand-derived pullback -- so
# agreement to machine precision means the reverse pass is differentiable in
# the same sense the primal is, which is the whole claim.
#
# It is also worth saying what this replaces. The R side previously took a
# central finite difference of the same gradient at step 1e-3; that agrees with
# the exact Hessian to roughly 1e-6 relative, and its error depends on a step
# that no single value can suit for a parameter vector mixing log standard
# deviations with unconstrained correlations. Both are asserted below, because
# "the new one is more accurate" is the reason for its existence.

function _hessian_test_model()
    df = DataFrame(
        matrix = [:T0MEANS, :T0MEANS, :LAMBDA, :LAMBDA, :LAMBDA, :LAMBDA,
                  :DRIFT, :DRIFT, :DRIFT, :DRIFT,
                  :DIFFUSION, :DIFFUSION, :DIFFUSION, :DIFFUSION,
                  :MANIFESTVAR, :MANIFESTVAR, :MANIFESTVAR, :MANIFESTVAR,
                  :MANIFESTMEANS, :MANIFESTMEANS, :CINT, :CINT,
                  :T0VAR, :T0VAR, :T0VAR, :T0VAR,
                  :JAx, :JAx, :JAx, :JAx, :Jy, :Jy, :Jy, :Jy, :PARS],
        row = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
               1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
        col = [1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1,
               1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1],
        parnumber = Union{Missing,Int}[
            1, 2, missing, missing, missing, missing,
            3, 4, 5, 6,
            7, 8, missing, 9,
            missing, missing, missing, missing,
            10, 11, 12, 13,
            missing, missing, missing, missing,
            missing, missing, missing, missing,
            missing, missing, missing, missing, missing],
        value = Union{Missing,Float64}[
            missing, missing, 1.0, 0.0, 0.0, 1.0,
            missing, missing, missing, missing,
            missing, missing, 0.0, missing,
            0.1, 0.0, 0.0, 0.1,
            missing, missing, missing, missing,
            1.0, 0.0, 0.0, 1.0,
            missing, missing, missing, missing,
            1.0, 0.0, 0.0, 1.0, 0.0],
        transform = Union{Missing,String}[
            "param[1]", "param[2]", missing, missing, missing, missing,
            "-log1p_exp(param[3])", "param[4]", "param[5]", "-log1p_exp(param[6])",
            "log1p_exp(param[7])", "param[8]", missing, "log1p_exp(param[9])",
            missing, missing, missing, missing,
            "param[10]", "param[11]", "param[12]", "param[13]",
            missing, missing, missing, missing,
            missing, missing, missing, missing,
            missing, missing, missing, missing, missing],
        # State-dependent DRIFT cells, so the tape's group machinery -- the part
        # whose replay scratch had to stop being Float64 for any of this to work
        # -- is actually exercised.
        predicttransform = Union{Missing,String}[
            fill(missing, 26)...,
            "DRIFT[1,1]", "DRIFT[2,1]", "DRIFT[1,2]", "DRIFT[2,2]",
            fill(missing, 5)...],
        updatetransform = fill(missing, 35),
    )
    sp = ekf_from_data_frame(df)

    nsubjects = 7
    subject_starts = Int[]
    times = Float64[]
    columns = Vector{Float64}[]
    position = 1
    for s in 1:nsubjects
        push!(subject_starts, position)
        nobs = 3 + (s % 3)
        for t in 1:nobs
            push!(times, 0.6 * (t - 1))
            push!(columns, [0.4 * sin(s + t), 0.3 * cos(2s - t)])
            position += 1
        end
    end
    objective = ctsem_objective(sp, subject_starts, times, reduce(hcat, columns))
    values = [0.2, -0.1, 0.3, 0.05, -0.05, 0.25, -0.2, 0.1, -0.15,
              0.02, -0.03, 0.04, -0.06]
    return objective, values
end

@testset "forward-over-reverse Hessian matches ForwardDiff on the primal" begin
    objective, values = _hessian_test_model()

    H = ctsem_hessian(objective, values)
    reference = ForwardDiff.hessian(objective, values)

    @test size(H) == (length(values), length(values))
    @test all(isfinite, H)
    @test H ≈ reference rtol = 1e-10
    # Exactly symmetric, not merely nearly so: the two triangles are averaged.
    @test H == transpose(H)
end

@testset "the forward-over-reverse Hessian beats the finite difference it replaces" begin
    objective, values = _hessian_test_model()
    reference = ForwardDiff.hessian(objective, values)

    step = 1e-3
    finite = zeros(length(values), length(values))
    for i in eachindex(values)
        plus = copy(values); minus = copy(values)
        plus[i] += step; minus[i] -= step
        finite[:, i] = (ctsem_adjoint_gradient(objective, plus).gradient .-
                        ctsem_adjoint_gradient(objective, minus).gradient) ./ (2step)
    end
    finite = (finite .+ transpose(finite)) ./ 2

    exact_error = norm(ctsem_hessian(objective, values) - reference) / norm(reference)
    finite_error = norm(finite - reference) / norm(reference)
    @test exact_error < 1e-12
    @test finite_error > 1e-9          # the difference is real, not noise
    @test exact_error < finite_error / 1e3
end

@testset "differentiating the gradient leaves the ordinary gradient intact" begin
    # The adjoint workspace is cached on the objective and keyed on the scalar
    # type, so a Hessian call replaces a Float64 workspace with a dual one. The
    # next ordinary gradient has to rebuild it and be unaffected -- otherwise
    # asking for uncertainty would quietly corrupt every later evaluation.
    objective, values = _hessian_test_model()

    before = ctsem_adjoint_gradient(objective, values)
    ctsem_hessian(objective, values)
    after = ctsem_adjoint_gradient(objective, values)

    @test after.value == before.value
    @test after.gradient == before.gradient
end

@testset "the Hessian is chunk-size independent" begin
    objective, values = _hessian_test_model()
    full = ctsem_hessian(objective, values; chunk=length(values))
    @test ctsem_hessian(objective, values; chunk=1) ≈ full rtol = 1e-12
    @test ctsem_hessian(objective, values; chunk=4) ≈ full rtol = 1e-12
end
