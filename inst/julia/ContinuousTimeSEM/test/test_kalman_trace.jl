# Per-row Kalman output.
#
# The point of the design is that prediction rides on the pass the likelihood
# already makes, so the first thing worth testing is that it still *is* that
# pass: the traced run must produce the same likelihood, row by row, and the
# adjoint must be unaffected by the new hooks existing. After that, the checks
# are the identities that define a filter and a smoother, rather than a
# comparison against the code that produced the numbers.

function _kalman_test_setup()
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
        predicttransform = Union{Missing,String}[
            fill(missing, 26)...,
            "DRIFT[1,1]", "DRIFT[2,1]", "DRIFT[1,2]", "DRIFT[2,2]",
            fill(missing, 5)...],
        updatetransform = fill(missing, 35),
    )
    sp = ekf_from_data_frame(df)
    times = [0.0, 0.5, 1.2, 2.0, 0.0, 0.4, 1.0]
    data = reduce(hcat, [[0.3, -0.1], [0.2, 0.4], [NaN, 0.1], [NaN, NaN],
        [0.5, 0.0], [0.1, -0.3], [0.2, 0.2]])
    objective = ctsem_objective(sp, [1, 5], times, data)
    values = [0.2, -0.1, 0.3, 0.05, -0.05, 0.25, -0.2, 0.1, -0.15,
              0.02, -0.03, 0.04, -0.06]
    return (sp=sp, objective=objective, values=values, data=data)
end

@testset "the traced pass is the same pass" begin
    setup = _kalman_test_setup()
    k = ctsem_kalman(setup.objective, setup.values)

    # Same likelihood, and the same likelihood row by row -- a trace that
    # perturbed the filter would show up here before anywhere else.
    plain = ctsem_evaluate(setup.objective, setup.values; gradient=false).value
    @test sum(k.subject_loglik) ≈ plain rtol = 1e-14
    @test sum(k.llrow) ≈ plain rtol = 1e-12

    # The new hooks are no-ops for the adjoint tape, so the gradient is
    # untouched by their existence.
    adjoint = ctsem_adjoint_gradient(setup.objective, setup.values)
    @test adjoint.value ≈ plain rtol = 1e-12
    forward = ctsem_evaluate(setup.objective, setup.values; gradient=true,
        gradient_method="forward")
    @test adjoint.gradient ≈ forward.gradient rtol = 1e-8

    @test size(k.eta) == (3, 7, 2)
    @test size(k.etacov) == (3, 7, 2, 2)
    @test size(k.y) == (3, 7, 2)
    @test size(k.ycov) == (3, 7, 2, 2)
    @test k.subject == [1, 1, 1, 1, 2, 2, 2]
end

@testset "filter identities hold at every row" begin
    setup = _kalman_test_setup()
    k = ctsem_kalman(setup.objective, setup.values)
    layout = ctsem_parameter_layout(setup.objective)
    flat = ctsem_parameter_matrices(setup.objective, setup.values)
    block(name) = begin
        j = findfirst(==(name), layout.matrix)
        reshape(flat[(layout.offset[j]+1):(layout.offset[j]+layout.nrow[j]*layout.ncol[j]), 1],
            layout.nrow[j], layout.ncol[j])
    end

    # A subject's first row starts from the population prior.
    for first in (1, 5)
        @test k.eta[1, first, :] ≈ vec(block("T0MEANS"))
        @test k.etacov[1, first, :, :] ≈ block("T0cov")
    end

    # The update can only remove uncertainty, and must remove some wherever
    # something was observed.
    for r in 1:7
        prior = k.etacov[1, r, :, :]
        upd = k.etacov[2, r, :, :]
        @test tr(upd) <= tr(prior) + 1e-12
    end
    observed = [any(!isnan, setup.data[:, r]) for r in 1:7]
    @test tr(k.etacov[2, 1, :, :]) < tr(k.etacov[1, 1, :, :])

    # Row 4 has nothing observed: the update is the identity there, and it
    # contributes nothing to the likelihood.
    @test !observed[4]
    @test k.eta[2, 4, :] ≈ k.eta[1, 4, :]
    @test k.etacov[2, 4, :, :] ≈ k.etacov[1, 4, :, :]
    @test k.llrow[4] == 0

    # The manifest quantities are the ones the model implies from the state,
    # over *all* manifest variables including the unobserved ones.
    LAMBDA = block("LAMBDA")
    MANIFESTMEANS = block("MANIFESTMEANS")
    MANIFESTcov = block("MANIFESTcov")
    Jy = block("Jy")
    for r in 1:7, kind in 1:3
        @test k.y[kind, r, :] ≈ MANIFESTMEANS[:, 1] .+ LAMBDA * k.eta[kind, r, :]
        @test k.ycov[kind, r, :, :] ≈ Jy * k.etacov[kind, r, :, :] * Jy' .+ MANIFESTcov
    end
end

@testset "the smoother runs backwards over each subject" begin
    setup = _kalman_test_setup()
    k = ctsem_kalman(setup.objective, setup.values)

    # At a subject's last row there is nothing left to condition on, so the
    # smoothed estimate is the filtered one.
    for last in (4, 7)
        @test k.eta[3, last, :] ≈ k.eta[2, last, :]
        @test k.etacov[3, last, :, :] ≈ k.etacov[2, last, :, :]
    end

    # Everywhere else smoothing uses no less data than filtering, so it cannot
    # be less certain, and where later rows carry information it must differ.
    for r in 1:7
        @test tr(k.etacov[3, r, :, :]) <= tr(k.etacov[2, r, :, :]) + 1e-12
    end
    for r in (1, 2, 5, 6)
        @test !isapprox(k.eta[3, r, :], k.eta[2, r, :]; rtol=1e-8)
    end
    @test tr(k.etacov[3, 1, :, :]) < tr(k.etacov[2, 1, :, :])

    # Row 4 observes nothing, so it adds nothing for row 3 to be smoothed with:
    # the correction term is the smoothed-minus-prior difference at row 4, and
    # with no update there that difference is exactly zero.
    @test k.eta[3, 4, :] ≈ k.eta[1, 4, :]
    @test k.eta[3, 3, :] ≈ k.eta[2, 3, :]

    # Smoothing must not leak across the subject boundary: subject 2's first row
    # cannot depend on subject 1's data. Refiltering subject 2 alone has to give
    # the identical answer.
    alone = ctsem_objective(setup.sp, [1], [0.0, 0.4, 1.0], setup.data[:, 5:7])
    k2 = ctsem_kalman(alone, setup.values)
    @test k2.eta[3, :, :] ≈ k.eta[3, 5:7, :] rtol = 1e-10
    @test k2.etacov[3, :, :, :] ≈ k.etacov[3, 5:7, :, :] rtol = 1e-10
end

@testset "subject matrices carry the smoothed initial state" begin
    setup = _kalman_test_setup()
    k = ctsem_kalman(setup.objective, setup.values)
    layout = ctsem_parameter_layout(setup.objective)
    population = ctsem_parameter_matrices(setup.objective, setup.values)

    @test size(k.subject_matrices) == (layout.size, 2)
    t0 = findfirst(==("T0MEANS"), layout.matrix)
    for (i, first) in enumerate((1, 5))
        @test k.subject_matrices[(layout.offset[t0]+1):(layout.offset[t0]+2), i] ≈
            k.eta[3, first, :]
    end

    # This model has no individually varying parameters, so every matrix other
    # than T0MEANS must be the population one -- an individual-differences
    # readout that moved a fixed parameter would be reporting noise.
    for j in eachindex(layout.matrix)
        layout.matrix[j] == "T0MEANS" && continue
        rows = (layout.offset[j]+1):(layout.offset[j]+layout.nrow[j]*layout.ncol[j])
        for i in 1:2
            @test k.subject_matrices[rows, i] ≈ population[rows, 1]
        end
    end
end
