"""
Progress reporting for long fits.

# Why this can exist at all

Output written by the engine reaches the R console *while a call is still
running* -- JuliaConnectoR forwards the subprocess's stdout as it arrives rather
than collecting it at the end. Measured: a Julia function printing once a second
for six seconds had its first line visible in R after one second and its last
before the call returned. So a fit can report on itself as it goes, which is the
only kind of progress report worth having.

# What is reported, and what is not

The rule throughout is that a progress line should answer "is this going to
finish, and is it going anywhere sensible" -- not "what is the engine doing".
So the optimiser reports the objective and the gradient norm, because a log
posterior that has stopped moving while the gradient is still large is a fit in
trouble; and the sampler reports step size, tree depth and divergences, because
those are what say whether the draws will be usable before the draws exist.

Updates are on a *time* cadence rather than an iteration one. An iteration
cadence prints thousands of lines on a fast model and nothing on a slow one,
which is exactly backwards: the user wants reassurance at a roughly constant
rate whatever the model costs.

# Threads

With several chains running at once only the first reports, and says so. Lines
from four threads interleave into something unreadable, and a lock to prevent
that would serialise the very thing being reported on. One chain's progress is
representative when the chains are doing the same work.
"""

using Printf

"""Cursor to column zero; see `_emit`."""
const CARRIAGE = Char(13)   # cursor to column zero
const NEWLINE = Char(10)

"""
    CTSEMProgress(enabled; every, label)

A rate-limited reporter. `every` is the minimum seconds between lines.
"""
mutable struct CTSEMProgress
    enabled::Bool
    every::Float64
    started::Float64
    last::Float64
    label::String
    lines::Int
    overwrite::Bool
    width::Int
    # Highest iteration count already shown. A fit can run Optim twice -- the
    # backtracking fallback continues from where Hager-Zhang stopped -- and the
    # second run's callback counts from one again, so the printed counter ran
    # *backwards*: "prior warm-up 2/10" and then "prior warm-up 1/10", which
    # reads as a fit that lost its place. The work did continue, so the count
    # is held rather than rewound.
    shown::Int
end

# `overwrite` defaults on: a fit reports for as long as it runs, and one line per
# report turns a five-minute sample into hundreds of lines of scrollback that
# say nothing the last of them does not. The caller turns it off when the output
# is going somewhere a carriage return means nothing -- a log file, a knitr
# chunk, a non-interactive session -- where the same updates should be rarer and
# on their own lines instead.
CTSEMProgress(enabled::Bool; every::Real=0.0, label::String="",
    overwrite::Bool=true) =
    CTSEMProgress(enabled, every > 0 ? Float64(every) : (overwrite ? 0.4 : 5.0),
        time(), 0.0, label, 0, overwrite, 0, 0)

"""Seconds since this reporter was created."""
_elapsed(p::CTSEMProgress) = time() - p.started

"""Should a line be printed now? True at most once per `every` seconds."""
function _due(p::CTSEMProgress, force::Bool=false)
    p.enabled || return false
    now = time()
    (force || now - p.last >= p.every) || return false
    p.last = now
    return true
end

"""`1m 04s`, `12.4s`, `2h 05m` -- whichever reads best at that magnitude."""
function _duration(seconds::Real)
    seconds = max(0.0, Float64(seconds))
    seconds < 60 && return @sprintf("%.1fs", seconds)
    seconds < 3600 && return @sprintf("%dm %02ds", floor(Int, seconds / 60),
        round(Int, seconds % 60))
    return @sprintf("%dh %02dm", floor(Int, seconds / 3600),
        round(Int, (seconds % 3600) / 60))
end

"""
    CTSEMConvergence(tolerance)

How far the optimiser has come toward its stopping rule, as a percentage.

# The formula is ctsem's, not a new one

`sgd.R` has printed a "Progress est." for years on the stan path, and this is
the same expression: a log-scale interpolation between the worst value seen,
the current one, and the tolerance. What it reports is the fraction of the
*orders of magnitude* between the starting point and the criterion that have
been covered -- not a fraction of the time, and not a fraction of the
iterations.

# Why the gradient and not the change in log posterior

`sgd.R` feeds the formula a log-posterior change because sgd *stops* on one.
This optimiser does not: `f_tol` is 0 by default, so Optim's `f` criterion is
off, and the fit ends when `maxabs(g) < g_tol` or when it runs out of
iterations. The change in log posterior therefore has no tolerance here to be
measured against, and inventing one would be inventing a criterion rather than
reporting against the one in force. The gradient is the criterion, so the
gradient is what is reported -- and it is the same number `g_tol` is compared
with, which is why the header can state the rule and the line can state the
distance to it in the same terms.

# Why the reported value never goes down

`sgd.R` does not clamp its output, and does not need to: the quantity it feeds
in is already smoothed -- a difference of running maxima over an
`nconvergeiter` window, non-negative by construction. Optim hands the callback
a raw per-iteration `g_norm`, and a line search that has to lengthen its step
raises it. Measured over 1184 iterations of eight fits (a well-conditioned one,
two- and three-latent, a nearly unidentified one, a saturating one, a binary
one, the state-augmented route and the Laplace route), the unclamped percentage
stepped backwards on 40--51% of iterations of the fits that grind -- by up to
9.9 points, and on 101 of the 200 lines a user would actually have seen on the
state-augmented fit. A percentage that goes backwards on half its updates is
not readable as progress, so the best reached so far is what is reported. The
raw gradient is on the same line and is where a fit going wrong is visible;
nothing is hidden by holding this number.

# What it is worth, measured

Against the fraction of iterations actually done, the estimate is one-sided:
across the five fits that ran to a stopping decision it never over-claimed by
more than 5.3 points, and it under-claimed by up to 53. The under-claiming is
concentrated in the last few iterations of a *healthy* fit, where L-BFGS
converges superlinearly and crosses five or more orders of magnitude in a
handful of steps -- so a number linear in `log|g|` must lag there, and there is
no fixing that without abandoning the criterion it is measured against.

That lag matters least where it is largest. A fit that converges superlinearly
is a fit that is nearly over, and it is on screen for a second; the fits that
show this number for minutes are the ones that converge linearly, where
`log|g|` really is linear in the iteration count and the estimate tracked the
iteration fraction to within about five points throughout. And on the fit that
never converged at all -- 1000 iterations to the cap -- it sat at 25--28% and
never promised completion, which is the one thing an iteration fraction cannot
do.

So it is an estimate that runs late and never runs early, which is the safe
direction: the failure it cannot have is telling someone a fit is nearly done
when it is not.

# The test the removed time estimate failed

`_progress_optimise` records what the old "time remaining" said at the end: 13.6
seconds, 0.27 seconds before the fit finished. The same question, asked of this
number on fits long enough to print more than one line -- the last *periodic*
line, not the forced closing one, since that is what a user's screen holds while
the fit ends:

    fit                        iterations   last line   said     outcome
    3 latent, 400 subjects         41       iter 39     99.1%    converged
    Laplace, 400 subjects          11       iter  8     89.0%    converged
    state augmented, 200 subj    1000       iter 965    45.8%    hit the cap
    nearly unidentified            46       iter  0      0.0%    saturated

The first three are the answer. The fourth is the limit of the whole mechanism
rather than of this number: that fit ran its 46 iterations inside a single 0.4s
cadence interval, so one line printed and the closing line followed it
immediately. Nothing printed on a time cadence can report on a fit that finishes
inside one interval.

# Where it is not reported

Not on a stage whose cap is the plan. The prior warm-up takes exactly its ten
iterations, so `7/10` is already exact there and an estimate would replace a
correct denominator with a guessed one.
"""
mutable struct CTSEMConvergence
    tolerance::Float64
    # Running worst, as in `sgd.R`'s `lpdiff1`. A new worst rebases the scale
    # rather than sending the percentage out of range.
    worst::Float64
    # Best reported so far; see above. NaN until the first estimate exists.
    best::Float64
end

CTSEMConvergence(tolerance::Real) = CTSEMConvergence(Float64(tolerance), 0.0, NaN)

"""
    _convergence_percent!(c, current)

Record this iteration's criterion value and return the percentage to report, or
`NaN` when no honest estimate exists -- no positive tolerance to aim at, or a
value that is not a finite non-negative number.
"""
function _convergence_percent!(c::CTSEMConvergence, current::Real)
    value = Float64(current)
    (c.tolerance > 0 && isfinite(value) && value >= 0) || return NaN
    # At or inside the criterion. Reached before the log below, which would be
    # asked for `log(0)` at an exactly zero gradient.
    if value <= c.tolerance
        c.best = 100.0
        return 100.0
    end
    c.worst = max(c.worst, value)
    # `value > tolerance` and `worst >= value`, so the span is positive.
    reached = clamp(100.0 * (1 - log(value / c.tolerance) /
        log(c.worst / c.tolerance)), 0.0, 100.0)
    c.best = isnan(c.best) ? reached : max(c.best, reached)
    return c.best
end

"""
    _progress_line(p, done, total, fields...)

One line: what fraction is done, how fast, how long is left, then whatever the
caller thinks matters.

**Only for work with a target it will actually reach** -- the sampler, which
takes exactly `nwarmup` and then exactly `ndraws` iterations. There the
extrapolation is sound, and measured so: a 200-draw chain reported "0.2s at
this rate" at draw 188 and finished 0.2s later.

An optimiser is the opposite case and uses `_progress_optimise` instead. The
note there records what extrapolating against `maxiter` actually produced.
"""
function _progress_line(p::CTSEMProgress, done::Integer, total::Integer,
    fields::AbstractString...)
    elapsed = _elapsed(p)
    rate = done <= 0 ? 0.0 : done / elapsed
    remaining = (rate <= 0 || total <= done) ? 0.0 : (total - done) / rate
    p.lines += 1
    parts = [@sprintf("%s %5d/%-5d", p.label, done, total),
             @sprintf("%5.1f/s", rate),
             @sprintf("%7s at this rate", _duration(remaining))]
    # Empty fields are dropped rather than joined: a caller with a field that
    # only sometimes applies passes "" for it, and joining that leaves a
    # separator with nothing after it.
    _emit(p, "  " * join(vcat(parts, filter(!isempty, collect(fields))), " | "))
    return nothing
end

"""
    _progress_optimise(p, done, cap, fields...; budget=false, percent=NaN)

The optimiser's line: iterations done, how far toward the stopping rule,
elapsed, then the caller's fields. Deliberately no time remaining, and usually
no denominator.

# Why there is no time remaining

An optimiser stops when the gradient is small enough, not when it has taken
`maxiter` iterations, so `(maxiter - done) / rate` extrapolates against a
number the fit never approaches. Measured on a 2-latent, 60-subject fit that
converged at iteration 307 of a cap of 1000, comparing each printed estimate
with the time the fit actually had left:

    iter    rate    predicted    realised    ratio
      17    41.0/s      24.0s       5.59s     4.3x
      99    48.1/s      18.7s       3.94s     4.7x
     209    51.0/s      15.5s       1.90s     8.1x
     297    51.8/s      13.6s       0.27s    51.0x

Every estimate was too long, the error grew as the fit approached its optimum,
and the last thing the user saw was "13.6s at this rate" 0.27s before the fit
finished. That is worse than useless: it is the number that makes someone kill
a fit that was nearly done. Elapsed time and an iteration count cannot be
wrong, so that is what this prints.

The rate is contaminated too, and separably. The prior warm-up has a cap of 10
that it always exhausts -- a correct denominator -- and still predicted "1m 14s
at this rate" for a stage that finished 0.1s later, because the first iteration
pays Julia's compilation and no rate measured across it describes the rest.

# Why the denominator usually is not shown

`maxiter` is a safety limit, not a target. Showing "297/1000" invites reading a
30% completion that means nothing, when the fit was in fact about to stop. So
the count stands alone until the cap is close enough to be a real risk -- past
half of it -- at which point the denominator is the story and appears.

`budget=true` is for a stage whose cap *is* the plan and is always reached, the
prior warm-up being the only one. There "7/10" is honest at every iteration and
reaching 10 is not a failure.

# What `percent` replaced

The iteration rate, which used to sit between the elapsed time and the caller's
fields. Two reasons, and the second is the deciding one.

A rate is worth printing when it is going to be multiplied by a count of work
remaining, and this line refuses to print a count of work remaining -- so the
rate was the residue of the estimate that was removed, and `_convergence_percent!`
answers the question it was standing in for directly. Elapsed and the iteration
count still give it to anyone who wants it.

And the line has to fit. A line that wraps cannot be overwritten in place: a
carriage return goes to the start of the last visual row, so the update leaves
the earlier rows behind as debris. The marginal line was already 81 characters
and the Laplace one, which also carries an inner-mode count, 102. Adding a
field without removing one puts those at 92 and 113; trading the rate for the
estimate leaves them at 83 and 98.
"""
function _progress_optimise(p::CTSEMProgress, done::Integer, cap::Integer,
    fields::AbstractString...; budget::Bool=false, percent::Real=NaN)
    elapsed = _elapsed(p)
    # Never rewind; see `shown`.
    done = p.shown = max(p.shown, Int(done))
    p.lines += 1
    counter = (budget || (cap > 0 && 2 * done >= cap)) ?
        @sprintf("%5d/%-5d", done, cap) : @sprintf("%5d iter ", done)
    parts = [@sprintf("%s %s", p.label, counter)]
    # "est." is not decoration. It is the same word `sgd.R` prints, and it is
    # the difference between a number that is checkable and one that is a
    # promise; see `CTSEMConvergence` for what the estimate is and is not.
    isnan(percent) || push!(parts, @sprintf("%3.0f%% est.", percent))
    push!(parts, @sprintf("%8s", _duration(elapsed)))
    _emit(p, "  " * join(vcat(parts, filter(!isempty, collect(fields))), " | "))
    return nothing
end

"""
    _progress_break(p)

End the current in-place line so something else can print on its own.

The engine's own `verbose` messages go to the same stdout as the progress line,
and a carriage-returned line has no newline on it -- so they landed *inside*
it: "|g| 2.70e+02ctsem_optimize: Hager-Zhang stopped after 10 iteration(s)".
Anything that prints while a fit is running calls this first.
"""
function _progress_break(p::CTSEMProgress)
    (p.enabled && p.overwrite && p.lines > 0) || return nothing
    print(NEWLINE)
    flush(stdout)
    p.width = 0
    # Not zero: the cursor is already at the start of a fresh line, so the next
    # update must not open with `_emit`'s leading newline and leave a blank one.
    p.lines = 1
    return nothing
end

"""
    _emit(p, text)

One update, in place or on its own line.

In place means a carriage return and no newline, padded to the width of the
longest line written so far -- without the padding a shorter update leaves the
tail of a longer one behind it, which reads as garbage rather than as progress.
Carriage returns survive the trip to the R console intact, which is what makes
this possible at all.
"""
function _emit(p::CTSEMProgress, text::AbstractString)
    if !p.overwrite
        println(text)
        flush(stdout)
        return nothing
    end
    # The first in-place update opens with a newline. A carriage return only
    # returns to the start of the current line, and that line may already hold
    # something written by R: the engine's progress goes to stdout while R's
    # messages go to stderr, the two are interleaved by the terminal, and the
    # result was "...afterwards skip it.  optimise 39/1000" on one line. One
    # newline costs nothing and guarantees the progress line starts on its own.
    if p.lines <= 1
        print(NEWLINE)
    end
    p.width = max(p.width, length(text))
    print(CARRIAGE, rpad(text, p.width))
    flush(stdout)
    return nothing
end

"""A closing line, with the total cost and whatever the caller wants said."""
function _progress_done(p::CTSEMProgress, fields::AbstractString...)
    p.enabled || return nothing
    text = "  " * join(vcat([@sprintf("%s done in %s", p.label,
        _duration(_elapsed(p)))], collect(fields)), " | ")
    if p.overwrite
        # Replaces the last in-place update, then ends the line so whatever
        # comes next starts cleanly.
        print(CARRIAGE, rpad(text, p.width), NEWLINE)
        flush(stdout)
        p.width = 0
    else
        println(text)
        flush(stdout)
    end
    return nothing
end

"""
A heading, printed once before the work starts.

The optimiser's counter has no denominator, so this is where the stopping rule
goes: one line saying what the fit is waiting for, rather than a fraction on
every line implying it is waiting for `maxiter`.

`p.lines` is bumped so the first update does not add `_emit`'s leading newline
on top of this line's own.
"""
function _progress_header(p::CTSEMProgress, text::AbstractString)
    p.enabled || return nothing
    println("  " * text)
    flush(stdout)
    p.lines = max(p.lines, 1)
    return nothing
end

# ---------------------------------------------------------------------------
# The trace, and the callback out to R
#
# Two consumers, two costs, so two mechanisms rather than one.
#
# The *trace* is what a convergence plot is drawn from, and it wants every
# iteration. Recording one is a push onto a vector in Julia -- free -- and the
# whole thing crosses the bridge once at the end, where a thousand iterations
# of three numbers is 24 kB and takes no measurable time.
#
# The *callback* is an R function invoked while the fit is still running, which
# is the only way a front end can show anything live. Measured through
# JuliaConnectoR it costs 0.5 ms per call: negligible occasionally, and 17% of a
# three-second fit if called every iteration. So it fires on the same time
# cadence as the printed line rather than on the iteration count -- what a
# watcher wants is a roughly constant update rate, not one update per unit of
# the engine's internal work.

"""
    CTSEMTrace(keys...)

Per-iteration record of whatever the caller thinks describes progress. Keys are
fixed at construction so the vectors stay type-stable.
"""
struct CTSEMTrace
    iteration::Vector{Int}
    values::Dict{Symbol,Vector{Float64}}
    keys::Vector{Symbol}
end

CTSEMTrace(keys::Symbol...) = CTSEMTrace(Int[],
    Dict(k => Float64[] for k in keys), collect(keys))

"""Append one iteration. Values are positional, in the key order given."""
function _record!(t::CTSEMTrace, iteration::Integer, values::Real...)
    length(values) == length(t.keys) ||
        throw(ArgumentError("trace expects $(length(t.keys)) values"))
    push!(t.iteration, Int(iteration))
    for (k, v) in zip(t.keys, values)
        push!(t.values[k], Float64(v))
    end
    return nothing
end

"""
Flat vectors rather than a matrix or a table: JuliaConnectoR moves a
`Vector{Float64}` without ceremony, and the R side rebuilds a data frame from
them in one step.
"""
function _trace_result(t::CTSEMTrace)
    out = Dict{Symbol,Any}(:iteration => t.iteration)
    for k in t.keys
        out[k] = t.values[k]
    end
    return (; (k => out[k] for k in vcat(:iteration, t.keys))...)
end

"""
    _invoke_callback(f, fields...)

Call R, and never let its failure take the fit with it. A callback is a
reporting convenience; an error inside one -- a plotting device that has gone
away, a typo in a user's function -- must not lose an optimisation that is
otherwise fine. It is reported once and then disabled.
"""
mutable struct CTSEMCallback
    f::Any
    alive::Bool
    every::Float64
    last::Float64
end
CTSEMCallback(f; every::Real=0.4) =
    CTSEMCallback(f, f !== nothing, Float64(every), 0.0)

"""
Its own clock, deliberately not the printed line's.

The two were shared at first and the callback then never fired: `_due` checks
`enabled`, which follows `verbose`, so a caller who passed a callback and left
`verbose = 0` -- the obvious combination for a front end that draws rather than
prints -- got nothing at all. A callback is a programmatic consumer and has no
reason to depend on whether anything is being printed.
"""
function _callback_due(cb::CTSEMCallback, force::Bool=false)
    cb.alive || return false
    now = time()
    (force || now - cb.last >= cb.every) || return false
    cb.last = now
    return true
end

function _invoke_callback(cb::CTSEMCallback, values...; force::Bool=false)
    _callback_due(cb, force) || return nothing
    try
        cb.f(values...)
    catch err
        cb.alive = false
        println("  progress callback failed and was disabled: ", err)
        flush(stdout)
    end
    return nothing
end
