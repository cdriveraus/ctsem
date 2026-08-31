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
        time(), 0.0, label, 0, overwrite, 0)

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
    _progress_line(p, done, total, fields...)

One line: what fraction is done, how fast, how long is left, then whatever the
caller thinks matters. Time remaining is extrapolated from the rate so far,
which is honest for a sampler (iterations cost about the same) and optimistic
for an optimiser (later iterations are usually cheaper), so it is labelled as an
estimate rather than presented as a fact.
"""
function _progress_line(p::CTSEMProgress, done::Integer, total::Integer,
    fields::AbstractString...)
    elapsed = _elapsed(p)
    rate = done <= 0 ? 0.0 : done / elapsed
    remaining = (rate <= 0 || total <= done) ? 0.0 : (total - done) / rate
    p.lines += 1
    parts = [@sprintf("%s %5d/%-5d", p.label, done, total),
             @sprintf("%5.1f/s", rate),
             # "at this rate", not "left". For an optimiser the two are very
             # different: the run below stopped four iterations after reporting
             # "15m 06s left", because its gradient was already 9.28e-09 and it
             # was about to converge -- the extrapolation assumes every one of
             # `maxiter` iterations will be taken, and convergence is precisely
             # the thing that stops that being true.
             @sprintf("%7s at this rate", _duration(remaining))]
    # Empty fields are dropped rather than joined: a caller with a field that
    # only sometimes applies passes "" for it, and joining that leaves a
    # separator with nothing after it.
    _emit(p, "  " * join(vcat(parts, filter(!isempty, collect(fields))), " | "))
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

"""A heading, printed once before the work starts."""
function _progress_header(enabled::Bool, text::AbstractString)
    enabled || return nothing
    println(text)
    flush(stdout)
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
