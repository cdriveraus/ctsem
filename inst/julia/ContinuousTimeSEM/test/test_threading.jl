# Single- vs multi-threaded agreement for the subject reduction.
#
# This is the acceptance gate `docs/src/adjoint-roadmap.md` lists as outstanding
# under "Acceptance Gates and Cleanup": the objective is a sum over subjects, and
# threading it changes only the order of that sum, so serial and threaded results
# must agree to floating-point reassociation and nothing more.
#
# The comparison is to 1e-12 relative rather than bitwise, deliberately. A
# threaded reduction sums per-chunk partials instead of subjects in order, so the
# last bits legitimately differ; demanding bitwise equality would either fail
# spuriously or force a summation order that defeats the point. A dropped or
# double-counted subject moves the answer far more than 1e-12, which is what this
# is actually guarding against.
#
# When Julia is started single-threaded (the default, and what `Pkg.test()` uses
# unless JULIA_NUM_THREADS says otherwise) the threaded path still runs -- with
# one chunk -- so these tests exercise the chunking arithmetic and the per-chunk
# workspace handling regardless. The genuinely concurrent comparison is skipped
# with an informative message rather than silently passing.

@testset "chunk ranges partition the subjects exactly" begin
    for n in 1:13, k in 1:6
        ranges = ContinuousTimeSEM._ctsem_chunk_ranges(n, k)
        @test sum(length, ranges) == n
        @test reduce(vcat, collect.(ranges)) == collect(1:n)
        @test length(ranges) == min(k, n)
        # near-equal: no chunk is more than one subject larger than another
        @test maximum(length, ranges) - minimum(length, ranges) <= 1
    end
end

@testset "chunk count respects the cap and the subject count" begin
    original = ContinuousTimeSEM.ctsem_max_chunks().max_chunks
    try
        ContinuousTimeSEM.ctsem_set_max_chunks!(1)
        @test ContinuousTimeSEM._ctsem_nchunks(10) == 1
        ContinuousTimeSEM.ctsem_set_max_chunks!(4)
        # never more chunks than subjects, and never more than the process has
        # threads -- a cap cannot conjure parallelism that does not exist
        @test ContinuousTimeSEM._ctsem_nchunks(2) <= 2
        @test ContinuousTimeSEM._ctsem_nchunks(100) <= Threads.nthreads()
        @test ContinuousTimeSEM._ctsem_nchunks(0) == 1
    finally
        ContinuousTimeSEM.ctsem_set_max_chunks!(original)
    end
end

@testset "serial and threaded results agree" begin
    # Eight subjects of unequal length, so the partition is uneven and a chunk
    # boundary falls inside the data rather than neatly between equal blocks.
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

    nsubjects = 8
    subject_starts = Int[]
    times = Float64[]
    columns = Vector{Float64}[]
    position = 1
    for s in 1:nsubjects
        push!(subject_starts, position)
        nobs = 3 + (s % 3)            # 3, 4 or 5 observations
        for t in 1:nobs
            push!(times, 0.6 * (t - 1))
            push!(columns, [0.4 * sin(s + t), 0.3 * cos(2s - t)])
            position += 1
        end
    end
    data = reduce(hcat, columns)
    objective = ctsem_objective(sp, subject_starts, times, data)
    values = [0.2, -0.1, 0.3, 0.05, -0.05, 0.25, -0.2, 0.1, -0.15,
              0.02, -0.03, 0.04, -0.06]

    original = ContinuousTimeSEM.ctsem_max_chunks().max_chunks
    try
        ContinuousTimeSEM.ctsem_set_max_chunks!(1)
        serial_value = objective(values)
        serial_adjoint = ctsem_adjoint_gradient(objective, values)
        serial_forward = ForwardDiff.gradient(objective, values)
        @test isfinite(serial_value)

        if Threads.nthreads() < 2
            @info string("Julia is running with one thread, so the threaded ",
                "reduction cannot be compared against the serial one here. ",
                "Start Julia with `-t auto` (or JULIA_NUM_THREADS) to exercise ",
                "it; the chunking arithmetic above is checked either way.")
        end

        for chunks in (2, 3, 8, 16)
            ContinuousTimeSEM.ctsem_set_max_chunks!(chunks)
            threaded_value = objective(values)
            threaded_adjoint = ctsem_adjoint_gradient(objective, values)
            threaded_forward = ForwardDiff.gradient(objective, values)
            @test threaded_value ≈ serial_value rtol = 1e-12
            @test threaded_adjoint.value ≈ serial_adjoint.value rtol = 1e-12
            @test threaded_adjoint.gradient ≈ serial_adjoint.gradient rtol = 1e-12
            @test threaded_forward ≈ serial_forward rtol = 1e-12
        end
    finally
        ContinuousTimeSEM.ctsem_set_max_chunks!(original)
    end
end

@testset "a failing subject surfaces from a threaded evaluation" begin
    # The threading-specific hazard is a *swallowed* failure: if one chunk's
    # task dies and the reduction still returns, the optimizer accepts a
    # gradient that is silently missing some subjects. `Threads.@sync` must
    # propagate the task's exception, and `ctsem_optimize`'s `fg!` must still be
    # able to catch it and reject the trial point.
    #
    # Worth recording, because two attempts to write this test as a
    # *non-finite* result failed: pathological parameter values in this engine
    # tend to **throw** rather than return NaN. `log` of a negative raises a
    # DomainError, and an infinite drift raises `InexactError: Int64(Inf)` from
    # the matrix exponential's scaling step. That is why `fg!` wraps the whole
    # evaluation in `try`, and it is also why the old "invalid innovation trials
    # poison the objective" test was deleted rather than rewritten: the
    # non-finite path is real but hard to reach from actual parameter values,
    # and manufacturing it by workspace surgery tests nothing.
    df = DataFrame(
        matrix = [:T0MEANS, :LAMBDA, :DRIFT, :DIFFUSION, :MANIFESTVAR,
                  :MANIFESTMEANS, :CINT, :T0VAR, :JAx, :Jy, :PARS],
        row = fill(1, 11), col = fill(1, 11),
        parnumber = Union{Missing,Int}[missing, missing, 1, missing, missing,
            missing, missing, missing, missing, missing, missing],
        value = Union{Missing,Float64}[0.0, 1.0, missing, 0.2, 0.1, 0.0, 0.0,
            1.0, missing, 1.0, 0.0],
        transform = Union{Missing,String}[missing, missing,
            "1 / param[1]", missing, missing, missing, missing, missing,
            missing, missing, missing],
        predicttransform = Union{Missing,String}[fill(missing, 8)...,
            "DRIFT[1,1]", missing, missing],
        updatetransform = fill(missing, 11),
    )
    sp = ekf_from_data_frame(df)
    objective = ctsem_objective(sp, [1, 4], collect(0.0:0.5:2.5),
        reshape(collect(0.1:0.1:0.6), 1, 6))

    # `fg!` catches by wrapping the whole evaluation; mirror that here rather
    # than asserting a specific exception type, which differs between the
    # serial path (the error itself) and the threaded one (wrapped by @sync).
    #
    # A non-finite result counts as failure as much as a thrown one does. The
    # matrix exponential used to throw `InexactError` on an infinite drift and
    # now returns NaN deliberately, so that every caller can treat a bad trial
    # point as an invalid point rather than an exception -- which is what this
    # test's own comment above describes wanting. Asserting only on the throw
    # left this failing three times over once that change landed.
    failed(values) = try
        !isfinite(ctsem_adjoint_gradient(objective, values).value)
    catch
        true
    end

    original = ContinuousTimeSEM.ctsem_max_chunks().max_chunks
    try
        ContinuousTimeSEM.ctsem_set_max_chunks!(1)
        # State the premise as an assertion: a finite point must still work, so
        # a test that passes because *everything* fails is not possible here.
        @test isfinite(ctsem_adjoint_gradient(objective, [-2.0]).value)
        @test failed([0.0])

        for chunks in (2, 4)
            ContinuousTimeSEM.ctsem_set_max_chunks!(chunks)
            @test isfinite(ctsem_adjoint_gradient(objective, [-2.0]).value)
            @test failed([0.0])
        end
    finally
        ContinuousTimeSEM.ctsem_set_max_chunks!(original)
    end
end
