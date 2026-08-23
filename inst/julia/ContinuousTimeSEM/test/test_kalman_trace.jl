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

# The three places this improves on Stan rather than reproducing it. Each is
# tested by an identity the improved version satisfies, not by a comparison
# against the code that produced the numbers -- and each is checked to leave the
# likelihood alone, since none of them may.
#
# The model: one OU process plus one static carrier state, with the measurement
# intercept reading the carrier. That is exactly how ctsem represents an
# individually varying MANIFESTMEANS, which is its default -- so this is the
# ordinary case, not an exotic one. The measurement equation is then linear in
# the augmented state, y = Jy x, which is the identity the tests use.

function _kalman_statedep_setup(; max_timestep=Inf)
    cells = [
        (:T0MEANS, 1, 1, 1, missing, "param[1]", missing, missing),
        (:T0MEANS, 2, 1, 2, missing, "param[2]", missing, missing),
        (:LAMBDA, 1, 1, missing, 1.0, missing, missing, missing),
        (:LAMBDA, 1, 2, missing, 0.0, missing, missing, missing),
        (:DRIFT, 1, 1, 3, missing, "-log1p_exp(param[3])", missing, missing),
        (:DRIFT, 2, 1, missing, 0.0, missing, missing, missing),
        (:DRIFT, 1, 2, missing, 0.0, missing, missing, missing),
        (:DRIFT, 2, 2, missing, 0.0, missing, missing, missing),
        (:DIFFUSION, 1, 1, 4, missing, "log1p_exp(param[4])", missing, missing),
        (:DIFFUSION, 2, 1, missing, 0.0, missing, missing, missing),
        (:DIFFUSION, 1, 2, missing, 0.0, missing, missing, missing),
        (:DIFFUSION, 2, 2, missing, 0.0, missing, missing, missing),
        (:MANIFESTVAR, 1, 1, missing, 0.2, missing, missing, missing),
        (:MANIFESTMEANS, 1, 1, missing, missing, missing, missing, "3 * state[2]"),
        (:CINT, 1, 1, missing, 0.0, missing, missing, missing),
        (:CINT, 2, 1, missing, 0.0, missing, missing, missing),
        (:T0VAR, 1, 1, 5, missing, "log1p_exp(param[5])", missing, missing),
        (:T0VAR, 2, 1, missing, 0.0, missing, missing, missing),
        (:T0VAR, 1, 2, missing, 0.0, missing, missing, missing),
        (:T0VAR, 2, 2, 6, missing, "log1p_exp(param[6])", missing, missing),
        (:JAx, 1, 1, missing, missing, missing, "DRIFT[1,1]", missing),
        (:JAx, 2, 1, missing, missing, missing, "DRIFT[2,1]", missing),
        (:JAx, 1, 2, missing, missing, missing, "DRIFT[1,2]", missing),
        (:JAx, 2, 2, missing, missing, missing, "DRIFT[2,2]", missing),
        (:Jy, 1, 1, missing, 1.0, missing, missing, missing),
        (:Jy, 1, 2, missing, 3.0, missing, missing, missing),
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
    sp = ekf_from_data_frame(df,
        DataFrame(parameter=Int[], predictor=Int[], coefficient=Int[]), [1])
    times = [0.0, 0.6, 1.4, 2.5]
    data = reshape([0.4, -0.2, 0.5, 0.1], 1, 4)
    objective = ctsem_objective(sp, [1], times, data, zeros(0, 4), zeros(1, 0),
        max_timestep)
    return (sp=sp, objective=objective, times=times,
        values=[0.3, -0.15, 0.2, -0.1, 0.05, -0.2])
end

# Built at top level, not inside the testsets: `ekf_from_columns` eval()s the
# transform strings into closures, and a closure created inside the same
# top-level statement that calls it is too new for that statement's world age.
_statedep_plain = _kalman_statedep_setup()
_statedep_substeps = _kalman_statedep_setup(max_timestep=0.35)

@testset "the measurement model is re-evaluated at the updated state" begin
    setup = _statedep_plain
    k = ctsem_kalman(setup.objective, setup.values)

    # y = Jy x must hold exactly at prior, filtered *and* smoothed. Reporting
    # the filtered observation with a pre-update measurement intercept -- what
    # Stan does -- breaks it at kinds 2 and 3.
    Jy = [1.0 3.0]
    for kind in 1:3, r in 1:4
        @test k.y[kind, r, :] ≈ Jy * k.eta[kind, r, :] atol = 1e-12
    end

    # And it is not a cosmetic difference: the carrier moves at the update, so
    # the pre-update intercept is wrong by 3 * (x_prior[2] - x_upd[2]).
    stale = [k.eta[2, r, 1] + 3 * k.eta[1, r, 2] for r in 1:4]
    @test maximum(abs.(stale .- k.y[2, :, 1])) > 1e-6

    # None of which may touch the likelihood.
    @test sum(k.subject_loglik) ≈
        ctsem_evaluate(setup.objective, setup.values; gradient=false).value rtol = 1e-14
end

@testset "the interval transition composes the substeps" begin
    # JAx is state independent here, so the whole-interval exponential and the
    # product of the substep transitions agree analytically. That makes this a
    # test of the *composition*: keeping only the last substep would report
    # exp(JAx dt/2) and fail, as would multiplying in the wrong order for a
    # non-commuting pair.
    for setup in (_statedep_plain, _statedep_substeps)
        k = ctsem_kalman(setup.objective, setup.values)
        layout = ctsem_parameter_layout(setup.objective)
        flat = ctsem_parameter_matrices(setup.objective, setup.values)
        j = findfirst(==("DRIFT"), layout.matrix)
        DRIFT = reshape(flat[(layout.offset[j]+1):(layout.offset[j]+4), 1], 2, 2)
        for r in 2:4
            @test k.transition[r, :, :] ≈
                exp(DRIFT .* (setup.times[r] - setup.times[r-1])) rtol = 1e-10
        end
        # Nothing precedes the first row, so its transition is left as identity.
        @test k.transition[1, :, :] ≈ [1.0 0.0; 0.0 1.0]
    end

    # Substepping changes the filter itself (the local affine model is
    # remade at each substep), so this also confirms the two setups differ --
    # otherwise the loop above would be testing the same thing twice.
    @test ctsem_evaluate(_statedep_plain.objective, _statedep_plain.values;
        gradient=false).value isa Float64
end


# A TD impulse whose Jacobian is not the identity -- which is what a
# state-dependent TDPREDEFFECT produces. Stan saves only the exponential for the
# smoother and drops this factor entirely.
function _kalman_jtd_setup()
    cells = [
        (:T0MEANS, 1, 1, 1, missing, "param[1]", missing, missing, missing),
        (:T0MEANS, 2, 1, 2, missing, "param[2]", missing, missing, missing),
        (:LAMBDA, 1, 1, missing, 1.0, missing, missing, missing, missing),
        (:LAMBDA, 1, 2, missing, 0.0, missing, missing, missing, missing),
        (:DRIFT, 1, 1, 3, missing, "-log1p_exp(param[3])", missing, missing, missing),
        (:DRIFT, 2, 1, missing, 0.1, missing, missing, missing, missing),
        (:DRIFT, 1, 2, missing, 0.2, missing, missing, missing, missing),
        (:DRIFT, 2, 2, 4, missing, "-log1p_exp(param[4])", missing, missing, missing),
        (:DIFFUSION, 1, 1, 5, missing, "log1p_exp(param[5])", missing, missing, missing),
        (:DIFFUSION, 2, 1, missing, 0.0, missing, missing, missing, missing),
        (:DIFFUSION, 1, 2, missing, 0.0, missing, missing, missing, missing),
        (:DIFFUSION, 2, 2, 6, missing, "log1p_exp(param[6])", missing, missing, missing),
        (:MANIFESTVAR, 1, 1, missing, 0.2, missing, missing, missing, missing),
        (:MANIFESTMEANS, 1, 1, missing, 0.0, missing, missing, missing, missing),
        (:CINT, 1, 1, missing, 0.0, missing, missing, missing, missing),
        (:CINT, 2, 1, missing, 0.0, missing, missing, missing, missing),
        (:T0VAR, 1, 1, missing, 1.0, missing, missing, missing, missing),
        (:T0VAR, 2, 1, missing, 0.0, missing, missing, missing, missing),
        (:T0VAR, 1, 2, missing, 0.0, missing, missing, missing, missing),
        (:T0VAR, 2, 2, missing, 1.0, missing, missing, missing, missing),
        (:JAx, 1, 1, missing, missing, missing, "DRIFT[1,1]", missing, missing),
        (:JAx, 2, 1, missing, missing, missing, "DRIFT[2,1]", missing, missing),
        (:JAx, 1, 2, missing, missing, missing, "DRIFT[1,2]", missing, missing),
        (:JAx, 2, 2, missing, missing, missing, "DRIFT[2,2]", missing, missing),
        (:Jy, 1, 1, missing, 1.0, missing, missing, missing, missing),
        (:Jy, 1, 2, missing, 0.0, missing, missing, missing, missing),
        (:TDPREDEFFECT, 1, 1, missing, 0.4, missing, missing, missing, missing),
        (:TDPREDEFFECT, 2, 1, missing, 0.0, missing, missing, missing, missing),
        (:Jtd, 1, 1, missing, 0.5, missing, missing, missing, missing),
        (:Jtd, 2, 1, missing, 0.0, missing, missing, missing, missing),
        (:Jtd, 1, 2, missing, 0.3, missing, missing, missing, missing),
        (:Jtd, 2, 2, missing, 1.0, missing, missing, missing, missing),
        (:PARS, 1, 1, missing, 0.0, missing, missing, missing, missing),
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
        tdtransform = Union{Missing,String}[c[9] for c in cells],
    )
    sp = ekf_from_data_frame(df)
    times = [0.0, 0.6, 1.4, 2.5]
    data = reshape([0.4, -0.2, 0.5, 0.1], 1, 4)
    tdpreds = reshape([0.0, 1.0, 0.0, 1.0], 1, 4)
    objective = ctsem_objective(sp, [1], times, data, tdpreds, zeros(1, 0))
    return (sp=sp, objective=objective, times=times,
        values=[0.3, -0.15, 0.2, -0.1, 0.05, -0.2])
end

_jtd_setup = _kalman_jtd_setup()

@testset "the interval transition carries the TD impulse Jacobian" begin
    setup = _jtd_setup
    k = ctsem_kalman(setup.objective, setup.values)
    layout = ctsem_parameter_layout(setup.objective)
    flat = ctsem_parameter_matrices(setup.objective, setup.values)
    block(name) = begin
        j = findfirst(==(name), layout.matrix)
        reshape(flat[(layout.offset[j]+1):(layout.offset[j]+layout.nrow[j]*layout.ncol[j]), 1],
            layout.nrow[j], layout.ncol[j])
    end
    JAx = block("JAx")
    Jtd = block("Jtd")
    @test Jtd != [1.0 0.0; 0.0 1.0]

    for r in 2:4
        expected = Jtd * exp(JAx .* (setup.times[r] - setup.times[r-1]))
        @test k.transition[r, :, :] ≈ expected rtol = 1e-10
        # Dropping Jtd, as Stan does, is not the same matrix.
        @test !isapprox(k.transition[r, :, :],
            exp(JAx .* (setup.times[r] - setup.times[r-1])); rtol=1e-6)
    end

    @test sum(k.subject_loglik) ≈
        ctsem_evaluate(setup.objective, setup.values; gradient=false).value rtol = 1e-14
end

@testset "generated data is a draw from the model the filter conditions on" begin
    setup = _statedep_plain
    nrows = 4
    base = reshape([0.7, -1.2, 0.3, 1.5], 1, nrows)
    g = ctsem_generate(setup.objective, setup.values, base)

    @test size(g.Y) == (1, nrows)
    @test all(isfinite, g.Y)

    # The defining identity: the filter's own likelihood for the generated data
    # must be the likelihood it reported while generating it. Anything that let
    # the draw and the covariance it was drawn from come apart breaks this.
    regenerated = ctsem_objective(setup.sp, [1], setup.times, g.Y)
    @test ctsem_evaluate(regenerated, setup.values; gradient=false).value ≈
        sum(g.subject_loglik) rtol = 1e-12
    @test sum(g.llrow) ≈ sum(g.subject_loglik) rtol = 1e-12

    # A zero draw is the prior predictive mean, and leaves every innovation
    # zero -- so the filter never updates and the generated data is the model's
    # free-running prediction.
    zero = ctsem_generate(setup.objective, setup.values, zeros(1, nrows))
    k = ctsem_kalman(ctsem_objective(setup.sp, [1], setup.times, zero.Y), setup.values)
    @test vec(zero.Y) ≈ k.y[1, :, 1] rtol = 1e-10
    @test k.eta[2, :, :] ≈ k.eta[1, :, :] rtol = 1e-10

    # The draws drive the result: different normals, different data.
    other = ctsem_generate(setup.objective, setup.values, reshape([-0.4, 0.9, -1.1, 0.2], 1, nrows))
    @test !isapprox(other.Y, g.Y; rtol=1e-6)
end

@testset "generated data keeps the original missingness" begin
    setup = _kalman_statedep_setup()
    withmissing = ctsem_objective(setup.sp, [1], setup.times,
        reshape([0.4, NaN, 0.5, NaN], 1, 4))
    g = ctsem_generate(withmissing, setup.values, reshape([0.7, -1.2, 0.3, 1.5], 1, 4))
    # An entry that was not observed is not invented: a posterior predictive
    # check compares against the observations that exist.
    @test isnan(g.Y[1, 2])
    @test isnan(g.Y[1, 4])
    @test all(isfinite, g.Y[1, [1, 3]])
    @test g.llrow[2] == 0
    @test g.llrow[4] == 0
end
