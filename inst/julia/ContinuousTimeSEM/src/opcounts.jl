"""
Operation counters for the discretisation kernels.

Deterministic instrumentation. Every matrix exponential, Lyapunov
factorisation and Fréchet block the engine forms increments one of these, and
the two `DiscretizationCache` guards record their hits and misses. They exist so
a test or a benchmark can ask "how many exponentials did that gradient cost?"
and get the same answer on every machine, which a wall-clock measurement on a
shared box cannot give. A `Ref{Int}` increment is a few instructions against
the O(n^3) kernel it sits beside, so they stay on.

`exp` counts every `my_exp!` call, which includes the Fréchet blocks that
`_ctsem_expm` routes through `my_exp!` (those at or below
`_CTSEM_SMALL_CHOLESKY`); `frechet` counts the blocks themselves, whichever
exponential they reach. `lyap_schur` counts factorisations actually performed,
not solves against cached factors. `series_intercept` and `series_noise` count
evaluations of the series kernels in `series_discretization.jl`, forward or
recomputed by a pullback alike; they stay at zero unless some interval was
short against a pivot of the drift. `series_tables` counts the series tables
built (`SeriesTable`): one per drift a run of intervals shares, so its ratio to
the kernels' counts says how far they were shared.

The `probe_*` counts are the Laplace inner solve's mode probes
(`_laplace_probe_modes`): units probed, probe starts, unit gradients spent
screening them, starts the screen found non-Gaussian along their line, Newton
solves from a start, and probes whose mode replaced the one found. On the
affect pilot one value probed every unit at about 100 starts each, nearly all
of them refused by the first screening gradient, and kept 2 better modes from
9 Newton solves: the screen's gradients were most of the evaluation.

A threaded subject loop increments these racily, so the counts are then lower
bounds. Read them from a single-threaded run when the exact number matters.
"""
const _CTSEM_OPCOUNT = (
    exp = Ref(0),
    exp_cache_hit = Ref(0),
    exp_cache_miss = Ref(0),
    lyap_schur = Ref(0),
    lyap_ksolve = Ref(0),
    lyap_cache_hit = Ref(0),
    lyap_cache_miss = Ref(0),
    frechet = Ref(0),
    series_intercept = Ref(0),
    series_noise = Ref(0),
    series_tables = Ref(0),
    probe_units = Ref(0),
    probe_starts = Ref(0),
    probe_gradients = Ref(0),
    probe_nongaussian = Ref(0),
    probe_newtons = Ref(0),
    probe_better = Ref(0),
)

"""Zero every operation counter."""
function ctsem_reset_opcounts!()
    for r in _CTSEM_OPCOUNT
        r[] = 0
    end
    return nothing
end

"""Current operation counts as a NamedTuple of integers."""
ctsem_opcounts() = map(r -> r[], _CTSEM_OPCOUNT)

export ctsem_opcounts, ctsem_reset_opcounts!
