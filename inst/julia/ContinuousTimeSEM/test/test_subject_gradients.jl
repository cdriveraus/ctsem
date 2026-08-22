# Per-subject gradient contributions -- the score matrix the R side's
# `scorecalc()` produces for the Stan backend, and what the OPG, sandwich and
# score-bootstrap uncertainty methods consume.
#
# The check that matters is that the rows sum to the summed gradient. The two
# routes share every primitive and differ only in where they accumulate, so a
# disagreement means one of the two per-subject shortcuts in the summed path --
# the shared parameter layer, or the batched Frechet contribution -- has leaked
# across a subject boundary. That is exactly the failure mode worth guarding:
# both shortcuts exist precisely because the summed gradient does not need
# per-subject correctness.

@testset "subject gradients sum to the summed gradient" begin
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

    original = ContinuousTimeSEM.ctsem_max_chunks().max_chunks
    try
        for chunks in (1, 2, 4)
            ContinuousTimeSEM.ctsem_set_max_chunks!(chunks)
            summed = ctsem_adjoint_gradient(objective, values)
            per_subject = ctsem_subject_gradients(objective, values)
            @test size(per_subject.scores) == (nsubjects, length(values))
            @test per_subject.value ≈ summed.value rtol = 1e-12
            @test vec(sum(per_subject.scores, dims=1)) ≈ summed.gradient rtol = 1e-10
            # Each row is a real contribution, not a copy of the total or zero.
            @test all(any(!iszero, per_subject.scores[i, :]) for i in 1:nsubjects)
        end
    finally
        ContinuousTimeSEM.ctsem_set_max_chunks!(original)
    end
end

@testset "subject gradients match one-subject objectives" begin
    # An independent route to the same numbers: build a separate objective
    # containing only subject i and take its ordinary gradient. This does not
    # share the per-subject bookkeeping with `ctsem_subject_gradients`, so it
    # catches an error in the row *indexing* that a sum check cannot.
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
    times = [0.0, 0.4, 1.1, 0.0, 0.5, 1.2, 2.0, 0.0, 0.7]
    data = reshape([0.1, -0.2, 0.3, 0.05, 0.4, -0.1, 0.2, -0.3, 0.15], 1, 9)
    starts = [1, 4, 8]
    values = [0.1, 0.3, -0.2, 0.05]

    objective = ctsem_objective(sp, starts, times, data)
    per_subject = ctsem_subject_gradients(objective, values)

    stops = vcat(starts[2:end] .- 1, length(times))
    for i in eachindex(starts)
        range = starts[i]:stops[i]
        alone = ctsem_objective(sp, [1], times[range] .- 0.0, data[:, range])
        @test ctsem_adjoint_gradient(alone, values).gradient ≈
            per_subject.scores[i, :] rtol = 1e-10
    end
end
