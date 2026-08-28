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
             @sprintf("%7s left", _duration(remaining))]
    _emit(p, "  " * join(vcat(parts, collect(fields)), " | "))
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
