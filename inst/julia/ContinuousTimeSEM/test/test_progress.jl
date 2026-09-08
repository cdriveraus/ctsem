# The progress line's formatting rules (progress.jl).
#
# These exist because the rules are about what is *not* printed, and nothing
# else in the suite would notice them coming back. The optimiser's line used to
# extrapolate a time remaining against `maxiter` -- a cap the fit never reaches
# -- and on a fit that converged at iteration 307 of 1000 the last thing on
# screen was "13.6s at this rate" 0.27s before it finished, a 51x overestimate.
# It also showed "297/1000", inviting a completion fraction that meant nothing.
#
# So: `_progress_optimise` prints elapsed and never an estimate, and carries a
# denominator only when the cap is close enough to be a real risk. `_progress_line`
# keeps its estimate, because the sampler's total is reached exactly.

using Printf

# Capture what a reporter writes. `overwrite=false` keeps each update on its own
# line, which is what makes the text testable at all.
#
# `stderr`, because that is where the engine writes: R's `stderr()` is where
# `message()` goes, so a front end that styles R's messages styles these too.
# See `_console()` in progress.jl.
#
# Through a temp file rather than an IOBuffer: `redirect_stderr` takes a real
# stream, and an IOBuffer raises `MethodError: no method matching
# (::Base.RedirectStdStream)(::IOBuffer)`.
function _capture_progress(f::Function)
    path, io = mktemp()
    try
        redirect_stderr(io) do
            f()
            flush(stderr)
        end
        close(io)
        return read(path, String)
    finally
        isopen(io) && close(io)
        rm(path; force=true)
    end
end

_optimise_reporter(; label="optimise") =
    ContinuousTimeSEM.CTSEMProgress(true; label=label, overwrite=false)

@testset "the optimiser's line never estimates a time remaining" begin
    text = _capture_progress() do
        p = _optimise_reporter()
        ContinuousTimeSEM._progress_optimise(p, 297, 1000, "logpost -1232.52")
    end
    # The whole point: no extrapolation against a cap the fit will not reach.
    @test !occursin("at this rate", text)
    @test !occursin("left", text)
    # Elapsed and the count are what replace it, and neither can be wrong.
    @test occursin("297", text)
    @test occursin("logpost -1232.52", text)
    # The iteration rate went when the convergence estimate arrived: a rate is
    # for multiplying by a count of work remaining, and this line prints no
    # count of work remaining. See `_progress_optimise`.
    @test !occursin("/s", text)
    # And nothing is printed in its place unless a percentage is passed.
    @test !occursin("est.", text)
end

@testset "the sampler's line keeps its estimate" begin
    text = _capture_progress() do
        p = ContinuousTimeSEM.CTSEMProgress(true; label="sampling", overwrite=false)
        ContinuousTimeSEM._progress_line(p, 188, 200, "logp -253.32")
    end
    # `ndraws` is reached exactly, so extrapolating to it is sound here.
    @test occursin("at this rate", text)
    @test occursin("188/200", text)
end

@testset "the denominator appears only when the cap is a real risk" begin
    # Below half the cap: a bare count, because "34/1000" reads as 3% done when
    # the fit is in fact about to converge.
    below = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(), 34, 1000)
    end
    @test occursin("34 iter", below)
    @test !occursin("34/1000", below)

    # At and past half, running out of iterations is the story.
    at_half = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(), 500, 1000)
    end
    @test occursin("500/1000", at_half)
    @test !occursin("500 iter", at_half)
end

@testset "a budget stage always shows its denominator" begin
    # The prior warm-up runs its cap of 10 every time, so "3/10" is honest at
    # every iteration even though 3 is well below half.
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(label="prior warm-up"),
            3, 10; budget=true)
    end
    @test occursin("3/10", text)
    @test !occursin("3 iter", text)
    @test !occursin("at this rate", text)
end

@testset "the printed counter does not rewind" begin
    # The backtracking fallback is a second Optim run whose callback counts from
    # one again, which printed "prior warm-up 2/10" and then "1/10" -- a fit that
    # looks like it lost its place. The work did continue, so the count is held.
    p = _optimise_reporter()
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(p, 36, 1000)
        ContinuousTimeSEM._progress_optimise(p, 1, 1000)
        ContinuousTimeSEM._progress_optimise(p, 40, 1000)
    end
    lines = filter(!isempty, split(text, '\n'))
    @test length(lines) == 3
    @test occursin("36 iter", lines[1])
    # Held at 36 rather than rewound to 1.
    @test occursin("36 iter", lines[2])
    @test occursin("40 iter", lines[3])
    @test p.shown == 40
end

@testset "_progress_break ends an in-place line" begin
    # The engine's own `verbose` messages share a stream with the progress line
    # -- deliberately, see `_console()` -- and a carriage-returned line has no
    # newline on it, so they landed inside it: "|g| 2.70e+02ctsem_optimize:
    # Hager-Zhang stopped...".
    p = ContinuousTimeSEM.CTSEMProgress(true; label="optimise", overwrite=true)
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(p, 10, 1000)
        ContinuousTimeSEM._progress_break(p)
        print(ContinuousTimeSEM._console(), "a message")
    end
    @test endswith(text, "\na message")
    # Not reset to zero: the cursor is already on a fresh line, so the next
    # update must not add a leading newline and leave a blank one.
    @test p.lines == 1

    # A disabled or non-overwriting reporter has no open line to close.
    quiet = ContinuousTimeSEM.CTSEMProgress(false; label="optimise")
    @test _capture_progress() do
        ContinuousTimeSEM._progress_break(quiet)
    end == ""
end

@testset "the convergence estimate is sgd.R's formula on this optimiser's rule" begin
    # `sgd.R:405`:
    #   100 * (1 - log(current/tol) / log(worst/tol))
    # with `worst` the running maximum. Same expression here, with the gradient
    # norm in place of the log-posterior change because the gradient is what
    # this optimiser stops on -- `f_tol` is 0 by default, so there is no
    # log-posterior tolerance to interpolate against.
    c = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    reference(cur, worst, tol) = 100 * (1 - log(cur / tol) / log(worst / tol))
    # The first value is the worst by definition, so the fit starts at zero.
    @test ContinuousTimeSEM._convergence_percent!(c, 361.8) == 0.0
    @test ContinuousTimeSEM._convergence_percent!(c, 0.6775) ≈
        reference(0.6775, 361.8, 1e-8)
    @test ContinuousTimeSEM._convergence_percent!(c, 1.423e-8) ≈
        reference(1.423e-8, 361.8, 1e-8)
end

@testset "a new worst rebases the scale rather than leaving the range" begin
    # `sgd.R` does the same with `lpdiff1`: a value worse than any seen becomes
    # the new baseline, so the reported number stays a fraction of a span that
    # actually contains the current value.
    c = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    ContinuousTimeSEM._convergence_percent!(c, 1.447)
    ContinuousTimeSEM._convergence_percent!(c, 0.3)
    @test c.worst == 1.447
    # Iteration 2 of the nearly-unidentified fit: the line search lengthens its
    # step and the gradient overshoots every value seen so far.
    ContinuousTimeSEM._convergence_percent!(c, 3.387)
    @test c.worst == 3.387
    # Still in range, and still a percentage.
    @test 0 <= c.best <= 100
end

@testset "the reported percentage never goes backwards" begin
    # Optim hands the callback a raw per-iteration gradient norm, unlike the
    # window-smoothed quantity `sgd.R` feeds the same formula. Measured over
    # eight fits, the unclamped number stepped backwards on up to 51% of
    # iterations, by as much as 9.9 points. The best reached is what is shown.
    c = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    # The gradient trace of the two-latent fit's first six iterations, which
    # rises twice.
    seen = Float64[1.551, 1.043, 1.384, 0.7369, 1.074, 0.3661]
    reported = [ContinuousTimeSEM._convergence_percent!(c, g) for g in seen]
    @test all(diff(reported) .>= 0)
    # Held, not recomputed: iteration 3's 1.384 is worse than iteration 2's
    # 1.043, so the raw formula would have fallen there.
    @test reported[3] == reported[2]
    # The raw gradient is on the same line, so a fit going backwards is still
    # visible; only the estimate is held.
    @test reported[end] > reported[1]
end

@testset "the estimate is 100% at the criterion and needs no log of zero" begin
    c = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    ContinuousTimeSEM._convergence_percent!(c, 12.0)
    @test ContinuousTimeSEM._convergence_percent!(c, 1e-8) == 100.0
    # A gradient that underflows to exactly zero is the normal end of a healthy
    # fit, and `log(0/tol)` would be -Inf. Reached before the logarithm.
    d = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    ContinuousTimeSEM._convergence_percent!(d, 12.0)
    @test ContinuousTimeSEM._convergence_percent!(d, 0.0) == 100.0
    # And past it, rather than beyond 100.
    e = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    ContinuousTimeSEM._convergence_percent!(e, 12.0)
    @test ContinuousTimeSEM._convergence_percent!(e, 1e-14) == 100.0
end

@testset "no estimate is reported when there is nothing to estimate against" begin
    # A NaN gradient is a fit in trouble, not a fit at 0%. `NaN` here means
    # "say nothing", and `_progress_optimise` then prints no field at all --
    # the failure this whole line exists to avoid is a plausible number.
    c = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    @test isnan(ContinuousTimeSEM._convergence_percent!(c, NaN))
    @test isnan(ContinuousTimeSEM._convergence_percent!(c, Inf))
    @test isnan(ContinuousTimeSEM._convergence_percent!(c, -1.0))
    # `g_tol = 0` asks for an exactly zero gradient, which no span can be
    # measured against.
    @test isnan(ContinuousTimeSEM._convergence_percent!(
        ContinuousTimeSEM.CTSEMConvergence(0.0), 1.0))
    # A bad iteration does not destroy an estimate already earned.
    d = ContinuousTimeSEM.CTSEMConvergence(1e-8)
    ContinuousTimeSEM._convergence_percent!(d, 100.0)
    good = ContinuousTimeSEM._convergence_percent!(d, 1e-3)
    @test isnan(ContinuousTimeSEM._convergence_percent!(d, NaN))
    @test ContinuousTimeSEM._convergence_percent!(d, 1e-3) == good
end

@testset "the percentage appears on the line, labelled as an estimate" begin
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(), 297, 1000,
            "logpost -1232.52"; percent=38.4)
    end
    @test occursin("38% est.", text)
    # Where a reader looks for "how far": immediately after the counter, ahead
    # of the elapsed time.
    @test findfirst("38% est.", text)[1] < findfirst("logpost", text)[1]
    # A NaN percentage prints nothing rather than "NaN%".
    absent = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(), 297, 1000,
            "logpost -1232.52"; percent=NaN)
    end
    @test !occursin("est.", absent)
    @test !occursin("NaN", absent)
end

@testset "a budget stage keeps its exact fraction and gets no estimate" begin
    # The prior warm-up runs its cap of 10 every time, so "3/10" is already
    # correct. Replacing a right denominator with an estimated one is the
    # mistake this whole file exists to prevent, so the caller passes NaN and
    # the line carries the fraction alone.
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(label="prior warm-up"),
            3, 10; budget=true, percent=NaN)
    end
    @test occursin("3/10", text)
    @test !occursin("est.", text)
end

@testset "the line does not grow wide enough to wrap" begin
    # A wrapped line cannot be overwritten in place: a carriage return goes to
    # the start of the last visual row and leaves the earlier ones as debris.
    # The widest optimiser line is the Laplace one, which carries an inner-mode
    # count as well.
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(_optimise_reporter(), 297, 1000,
            @sprintf("logpost %11.2f", -1232.52),
            @sprintf("|g| %9.2e", 3.98e-13),
            @sprintf("inner %d/%d", 120, 120); percent=38.4)
    end
    @test maximum(length, split(strip(text, '\n'), '\n')) <= 100
end

@testset "durations read at the magnitude they are" begin
    @test ContinuousTimeSEM._duration(12.4) == "12.4s"
    @test ContinuousTimeSEM._duration(64) == "1m 04s"
    @test ContinuousTimeSEM._duration(7500) == "2h 05m"
    # Never negative, whatever the arithmetic upstream produced.
    @test ContinuousTimeSEM._duration(-5) == "0.0s"
end
