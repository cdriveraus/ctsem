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
# Through a temp file rather than an IOBuffer: `redirect_stdout` takes a real
# stream, and an IOBuffer raises `MethodError: no method matching
# (::Base.RedirectStdStream)(::IOBuffer)`.
function _capture_progress(f::Function)
    path, io = mktemp()
    try
        redirect_stdout(io) do
            f()
            flush(stdout)
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
    @test occursin("/s", text)
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
    # The engine's own `verbose` messages share stdout with the progress line,
    # and a carriage-returned line has no newline on it -- so they landed inside
    # it: "|g| 2.70e+02ctsem_optimize: Hager-Zhang stopped...".
    p = ContinuousTimeSEM.CTSEMProgress(true; label="optimise", overwrite=true)
    text = _capture_progress() do
        ContinuousTimeSEM._progress_optimise(p, 10, 1000)
        ContinuousTimeSEM._progress_break(p)
        print("a message")
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

@testset "the header states the stopping rule once" begin
    p = _optimise_reporter()
    text = _capture_progress() do
        ContinuousTimeSEM._progress_header(p,
            @sprintf("optimise: stops when |g| < %.0e, or at %d iterations", 1e-8, 1000))
    end
    @test occursin("stops when |g| < 1e-08", text)
    @test occursin("1000 iterations", text)
    # Bumped so the first update does not add `_emit`'s leading newline on top
    # of this line's own.
    @test p.lines >= 1
    # Nothing is printed when reporting is off.
    off = ContinuousTimeSEM.CTSEMProgress(false)
    @test _capture_progress() do
        ContinuousTimeSEM._progress_header(off, "anything")
    end == ""
end

@testset "durations read at the magnitude they are" begin
    @test ContinuousTimeSEM._duration(12.4) == "12.4s"
    @test ContinuousTimeSEM._duration(64) == "1m 04s"
    @test ContinuousTimeSEM._duration(7500) == "2h 05m"
    # Never negative, whatever the arithmetic upstream produced.
    @test ContinuousTimeSEM._duration(-5) == "0.0s"
end
