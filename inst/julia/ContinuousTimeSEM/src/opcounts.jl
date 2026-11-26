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
evaluations of the series kernels in `series_discretization.jl`, by the forward
pass (the adjoint's pullbacks read the chain it taped); they stay at zero unless some interval was
short against a pivot of the drift. `series_tables` counts the series tables
built (`SeriesTable`): one per drift a run of intervals shares, so its ratio to
the kernels' counts says how far they were shared.

The `probe_*` counts are the Laplace inner solve's mode probes
(`_laplace_probe_modes`): units probed, probe starts, unit values and unit
gradients spent screening them, starts the screen found non-Gaussian along
their line, Newton solves from a start, and probes whose mode replaced the one
found. On the affect pilot one value probed every unit at about 100 starts
each, nearly all of them refused by the first screening gradient, and kept 2
better modes from 9 Newton solves: the screen's gradients were most of the
evaluation.

The `inner_*` counts are the Laplace inner solve's Newton
(`_laplace_newton_unit_mode`): solves, their iterations, and the exact unit
curvatures they formed. Their ratios say how far a solve's start was from its
mode: an anchored solve one step from its anchor takes one iteration and one
curvature, a search from the origin several.

The `sweep_*` counts are the Laplace route's second-order member passes, a
member each: a unit curvature's (`_laplace_columns_jacobian`, a pass per lane
chunk), the seeded first-order sweeps (`_tan_laplace_seeded`) and the jets
(`_laplace_unit_jet_gradient`). A sweep costs about its member's pass whatever
its lanes, so these say where a unit's value and gradient go.

The `overshoot_*` counts are the overshoot probe's (`_ctsem_overshot`): trial
points, sets it skipped as already tried, full evaluations of a trial, the
screen's gradients at the estimate (`CTSEMProbeScreen`), its sampled
evaluations and the units in them, and trials it settled without a full
evaluation. A full evaluation costs every unit, so `overshoot_values * N +
overshoot_units` is the trials' cost in unit evaluations. The probe runs them
one after another, so these are exact at any thread count.

`fold_passes` is the member passes the folded initial means' mode solves take
(`_hybrid_fold_score!`, called by `_hybrid_fold_member_modes`): a pass is a
member's value and gradient at fixed deviations. It is atomic, so exact at any
thread count.

A threaded subject loop increments the others racily, so the counts are then
lower bounds. Read them from a single-threaded run when the exact number
matters.
"""
const _CTSEM_OPCOUNT = (
    exp = Ref(0),
    exp_cache_hit = Ref(0),
    exp_cache_miss = Ref(0),
    lyap_schur = Ref(0),
    lyap_cache_hit = Ref(0),
    lyap_cache_miss = Ref(0),
    frechet = Ref(0),
    series_intercept = Ref(0),
    series_noise = Ref(0),
    series_tables = Ref(0),
    probe_units = Ref(0),
    probe_starts = Ref(0),
    probe_values = Ref(0),
    probe_gradients = Ref(0),
    probe_nongaussian = Ref(0),
    probe_newtons = Ref(0),
    probe_better = Ref(0),
    inner_solves = Ref(0),
    inner_iterations = Ref(0),
    inner_curvatures = Ref(0),
    sweep_curvature = Ref(0),
    sweep_seeded = Ref(0),
    sweep_jet = Ref(0),
    overshoot_trials = Ref(0),
    overshoot_repeats = Ref(0),
    overshoot_values = Ref(0),
    overshoot_gradients = Ref(0),
    overshoot_samples = Ref(0),
    overshoot_units = Ref(0),
    overshoot_screened = Ref(0),
    fold_passes = Threads.Atomic{Int}(0),
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
