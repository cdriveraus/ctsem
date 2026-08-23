# Model-implied parameter matrices.
#
# These are what the R side's summary and plot functions are built on, so the
# checks worth making are the ones that would let a wrong summary look right:
#
#   * the flat layout and the values agree -- a block written at the wrong
#     offset would still produce a plausible-looking matrix of the right shape;
#   * the derived matrices satisfy their defining equations rather than merely
#     being reproduced by the same code that computed them;
#   * a state-dependent cell is reported at the state it was evaluated at, and
#     is named as state dependent, since a caller that treats it as a constant
#     is the failure mode this interface exists to prevent;
#   * a batch of parameter vectors gives column by column what one at a time
#     gives, because R always calls it in batch.

@testset "parameter matrices: layout, values and derived quantities" begin
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
    times = [0.0, 0.5, 1.2, 0.0, 0.4]
    data = reduce(hcat, [[0.3, -0.1], [0.2, 0.4], [-0.2, 0.1], [0.5, 0.0], [0.1, -0.3]])
    objective = ctsem_objective(sp, [1, 4], times, data)
    values = [0.2, -0.1, 0.3, 0.05, -0.05, 0.25, -0.2, 0.1, -0.15,
              0.02, -0.03, 0.04, -0.06]

    layout = ctsem_parameter_layout(objective)
    flat = ctsem_parameter_matrices(objective, values)
    @test size(flat) == (layout.size, 1)
    @test layout.size == sum(layout.nrow .* layout.ncol)
    @test layout.offset == vcat(0, cumsum(layout.nrow .* layout.ncol)[1:end-1])
    for name in ("DIFFUSIONcov", "MANIFESTcov", "T0cov", "asymDIFFUSIONcov", "asymCINT")
        @test name in layout.matrix
    end

    block(name) = begin
        k = findfirst(==(name), layout.matrix)
        reshape(flat[(layout.offset[k]+1):(layout.offset[k]+layout.nrow[k]*layout.ncol[k]), 1],
            layout.nrow[k], layout.ncol[k])
    end

    # The transformed values are the transforms, evaluated. Checking two of the
    # nonlinear ones is enough to catch a block written at the wrong offset,
    # which is the error a shape check cannot see.
    DRIFT = block("DRIFT")
    @test DRIFT[1, 1] ≈ -ContinuousTimeSEM.log1p_exp(values[3])
    @test DRIFT[2, 1] ≈ values[4]
    @test block("DIFFUSION")[1, 1] ≈ ContinuousTimeSEM.log1p_exp(values[7])
    @test block("T0MEANS") ≈ reshape(values[1:2], 2, 1)

    # Derived matrices, against their defining equations rather than against a
    # second call to the code that produced them.
    DIFFUSIONcov = block("DIFFUSIONcov")
    @test DIFFUSIONcov ≈ Matrix(ContinuousTimeSEM.sdcovsqrt2cov(block("DIFFUSION"), 0))
    @test issymmetric(round.(DIFFUSIONcov, digits=12))
    @test block("T0cov") ≈ Matrix(ContinuousTimeSEM.sdcovsqrt2cov(block("T0VAR"), 0))
    @test block("MANIFESTcov") ≈ Matrix(ContinuousTimeSEM.sdcovsqrt2cov(block("MANIFESTVAR"), 0))

    asym = block("asymDIFFUSIONcov")
    @test maximum(abs, DRIFT * asym + asym * DRIFT' + DIFFUSIONcov) < 1e-10
    @test maximum(abs, DRIFT * block("asymCINT") + block("CINT")) < 1e-10
end

@testset "state-dependent cells are evaluated at the given state and named" begin
    # JAx[1,1] follows the state rather than DRIFT, so the reported value has to
    # move with `state` -- and has to be flagged, because a caller treating it
    # as a constant is exactly what this interface exists to prevent.
    df = DataFrame(
        matrix = [:T0MEANS, :LAMBDA, :DRIFT, :DIFFUSION, :MANIFESTVAR,
                  :MANIFESTMEANS, :CINT, :T0VAR, :JAx, :Jy, :PARS],
        row = fill(1, 11), col = fill(1, 11),
        parnumber = Union{Missing,Int}[1, missing, 2, 3, missing, 4, missing,
            missing, missing, missing, missing],
        value = Union{Missing,Float64}[missing, 1.0, missing, missing, 0.1,
            missing, 0.0, 1.0, missing, 1.0, 0.0],
        transform = Union{Missing,String}["param[1]", missing,
            "-log1p_exp(param[2])", "log1p_exp(param[3])", missing, "param[4]",
            missing, missing, missing, missing, missing],
        predicttransform = Union{Missing,String}[fill(missing, 8)...,
            "DRIFT[1,1] - 0.5 * state[1]", missing, missing],
        updatetransform = fill(missing, 11),
    )
    sp = ekf_from_data_frame(df)
    objective = ctsem_objective(sp, [1], [0.0, 0.4, 1.1], reshape([0.1, -0.2, 0.3], 1, 3))
    values = [0.7, 0.3, -0.2, 0.05]

    layout = ctsem_parameter_layout(objective)
    @test layout.statedep_matrix == ["JAx"]
    @test layout.statedep_row == [1]
    @test layout.statedep_col == [1]

    jax = findfirst(==("JAx"), layout.matrix)
    drift = findfirst(==("DRIFT"), layout.matrix)
    at_default = ctsem_parameter_matrices(objective, values)
    at_state = ctsem_parameter_matrices(objective, values; state=[2.0])
    # The default state is T0MEANS, which here is values[1].
    @test at_default[layout.offset[jax]+1, 1] ≈
        at_default[layout.offset[drift]+1, 1] - 0.5 * values[1]
    @test at_state[layout.offset[jax]+1, 1] ≈
        at_default[layout.offset[drift]+1, 1] - 0.5 * 2.0
    # Only the state-dependent cell moved.
    moved = findall(k -> !isapprox(at_default[k, 1], at_state[k, 1]), 1:layout.size)
    @test moved == [layout.offset[jax] + 1]
end

@testset "a batch of parameter vectors matches one at a time" begin
    df = DataFrame(
        matrix = [:T0MEANS, :LAMBDA, :DRIFT, :DIFFUSION, :MANIFESTVAR,
                  :MANIFESTMEANS, :CINT, :T0VAR, :JAx, :Jy, :PARS],
        row = fill(1, 11), col = fill(1, 11),
        parnumber = Union{Missing,Int}[1, missing, 2, 3, missing, 4, missing,
            missing, missing, missing, missing],
        value = Union{Missing,Float64}[missing, 1.0, missing, missing, 0.1,
            missing, 0.0, 1.0, missing, 1.0, 0.0],
        transform = Union{Missing,String}["param[1]", missing,
            "-log1p_exp(param[2])", "log1p_exp(param[3])", missing, "param[4]",
            missing, missing, missing, missing, missing],
        predicttransform = Union{Missing,String}[fill(missing, 8)...,
            "DRIFT[1,1]", missing, missing],
        updatetransform = fill(missing, 11),
    )
    sp = ekf_from_data_frame(df)
    objective = ctsem_objective(sp, [1], [0.0, 0.4, 1.1], reshape([0.1, -0.2, 0.3], 1, 3))

    draws = [0.1 0.4 -0.3; 0.3 -0.2 0.5; -0.2 0.1 0.2; 0.05 -0.4 0.0]
    batched = ctsem_parameter_matrices(objective, draws)
    @test size(batched, 2) == 3
    for column in 1:3
        @test batched[:, column] ≈ ctsem_parameter_matrices(objective, draws[:, column])[:, 1]
    end
end

@testset "a non-stationary drift gives NaN asymptotics, not an error" begin
    # A positive drift has no asymptotic form. One such posterior draw must not
    # abort the summary of the others, so it comes back as NaN.
    df = DataFrame(
        matrix = [:T0MEANS, :LAMBDA, :DRIFT, :DIFFUSION, :MANIFESTVAR,
                  :MANIFESTMEANS, :CINT, :T0VAR, :JAx, :Jy, :PARS],
        row = fill(1, 11), col = fill(1, 11),
        parnumber = Union{Missing,Int}[1, missing, 2, 3, missing, 4, missing,
            missing, missing, missing, missing],
        value = Union{Missing,Float64}[missing, 1.0, missing, missing, 0.1,
            missing, 0.0, 1.0, missing, 1.0, 0.0],
        transform = Union{Missing,String}["param[1]", missing, "param[2]",
            "log1p_exp(param[3])", missing, "param[4]", missing, missing,
            missing, missing, missing],
        predicttransform = Union{Missing,String}[fill(missing, 8)...,
            "DRIFT[1,1]", missing, missing],
        updatetransform = fill(missing, 11),
    )
    sp = ekf_from_data_frame(df)
    objective = ctsem_objective(sp, [1], [0.0, 0.4, 1.1], reshape([0.1, -0.2, 0.3], 1, 3))
    layout = ctsem_parameter_layout(objective)
    k = findfirst(==("asymDIFFUSIONcov"), layout.matrix)

    stable = ctsem_parameter_matrices(objective, [0.1, -0.8, -0.2, 0.05])
    @test isfinite(stable[layout.offset[k]+1, 1])
    # DRIFT = 0 exactly: the Lyapunov system is singular.
    singular = ctsem_parameter_matrices(objective, [0.1, 0.0, -0.2, 0.05])
    @test isnan(singular[layout.offset[k]+1, 1])
end
