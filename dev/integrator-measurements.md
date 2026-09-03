# Integration and speed on the julia engine: measurements and decisions

Branch `integrator`, 2026-09-03. Every number below is from **dev1** (23-core
EPYC, Julia 1.12.7, `julia -t 1`, `BLAS.set_num_threads(1)`), taken while one
to two other jobs were running (load average 3 to 5). Times are the minimum of
seven repeats and are indicative; operation counts come from
`ctsem_opcounts()` and are exact. Scripts: `measure.jl`, `alloc_profile.jl`,
`lyap_threshold.jl` in the session scratchpad; the model builders are the same
`ekf_from_columns` path R uses.

## What the integrator is

Per prediction substep the filter forms `A = exp(JAx h)`, the intercept
`h φ1(JAx h) f(x)` through a linear solve, and the process covariance
`Q_d = X - A X A'` from the asymptotic Lyapunov solution `X`. The mean update is
exponential Rosenbrock-Euler exactly, so it is second order in the change of
`JAx` across the step; the covariance is the linearised recursion and is exact
for a linear model at any step. Substeps re-linearise at the step start. The
number of substeps is `ceil(dt / maxtimestep)` with one global `maxtimestep`,
default one step.

## Baseline, 200 subjects x 20 rows (100 x 20 for state-dependent)

`frechet` is the number of matrix-exponential Fréchet evaluations in one
adjoint gradient; `exp` the number of forward exponentials in one primal
evaluation. There are 3800 transitions (1900 with four substeps each for the
state-dependent models, so 7600 substeps).

| model | design | exp | frechet | adjoint ms | largest profile shares |
|---|---|---|---|---|---|
| linear n=2 | balanced | 0 | 1 | 22.7 | gc 30, fwd_update 17, rev_update 15 |
| linear n=2 | shared irregular (19 distinct dt) | 3800 | 3800 | 30.1 | frechet 16, rev_update 14 |
| linear n=2 | few distinct (3 values) | 2523 | 2536 | 27.7 | rev_update 16, gc 13, frechet 13 |
| linear n=2 | fully irregular | 3800 | 3800 | 30.2 | gc 15, frechet 15 |
| linear n=6 | balanced | 0 | 1 | 77.0 | rev_update 29, rev_predict 20 |
| linear n=6 | shared irregular | 3800 | 3800 | 105.5 | rev_update 22, frechet 16 |
| augmented CINT k=3 m=3 | shared irregular | 3800 | 3800 | 88.3 | frechet 19, rev_update 18 |
| augmented DRIFT k=2 m=4 | balanced | 3800 | 3800 | 99.1 | rev_group 20, frechet 15, gc 15 |
| state-dependent n=2, maxdt .25 | balanced | 7600 | 7600 | 49.5 | gc 40, frechet 14 |
| state-dependent n=6, maxdt .25 | balanced | 7600 | 7600 | 192.5 | lyap 22, rev_predict 18, frechet 17 |

Three things the counts say that the timings alone would not:

1. **Both caches were last-value.** The exp cache in `DiscretizationCache`
   held one `(JAx, dt)`; the Fréchet batch compared against the last pending
   direction. A shared wave schedule with only 19 distinct intervals therefore
   missed on every one of 3800 rows, and three interval lengths in random
   order missed two thirds of the time. Bundled data are mostly of this kind:
   AnomAuth has 2 distinct intervals, the ctExample sets 1 to 6, Oscillating
   53 of 2000; only ctstantestdat is fully irregular (270 of 270).
2. **Augmented models with individually varying drift never batch.** The
   predict group writes JAx cells, so the batch must flush before every group
   and `frechet = 3800` on a balanced panel. Individually varying CINT keeps
   JAx fixed and caches like a linear model.
3. **State-dependent models above four diffusing states paid a LAPACK Schur
   per substep**, twice in the adjoint (forward sweep and pullback), with 1 to
   5 KB allocated per call: 94% of the primal's allocation in that cell.

## Changes on this branch

- `my_exp_frechet!` (`frechet_exponential.jl`): the Al-Mohy and Higham (2009)
  recurrence computes `exp(A)` and `L(A, E)` together with about nineteen
  `n`-sized products, replacing the `2n x 2n` block exponential (six products
  at eight times the arithmetic each). Exact: same Padé approximant. The
  adjoint workspace owns the buffer, so a flush allocates nothing. Pinned
  against the block form at n = 1 to 20 and three norm scales
  (`test_frechet_exponential.jl`), and by the 96-assertion adjoint gate.
- Fréchet directions are batched in a table of up to 32 distinct
  `(JAx dt)` keys (`_frechet_enqueue!`), not only against the last one.
- The forward exponential cache is an `ExpTable` of distinct `dt` per `JAx`,
  and the reverse pass installs one table per chunk into every subject
  workspace it runs, so subjects share exponentials.
- The Lyapunov Schur threshold moved from 4 to 10 diffusing states
  (`_CTSEM_LYAP_SCHUR_ABOVE`), from a micro-benchmark on dev1 (packed solve
  faster to k = 10, crossover 11), and the packed buffer now caches its LU
  across calls as the Schur buffer already did. Without that cache the move
  made the linear 6-latent reverse pass 8% slower.
- Four per-row allocation sites removed: two runtime `Val` dispatches in the
  discrete-intercept solve, a runtime `Val` matvec in the measurement update,
  a per-row `zeros` on the `:theta` tape entry.
- `ctsem_opcounts()`: deterministic counters for all of the above.

## After (same cells, same machine, load 4 to 6 during this run)

Two runs: after the kernel, the two tables and the Schur threshold; and after
the allocation work (LU cache in the packed buffer, scratch for the two
pullbacks, the four `Val`/`zeros` sites). Counts are from the second run.

| model | design | exp | frechet | baseline ms | after tables | after allocation | alloc MB before / after |
|---|---|---|---|---|---|---|---|
| linear n=2 | balanced | 0 | 1 | 22.7 | 23.1 | 19.6 | 11.9 / 1.7 |
| linear n=2 | shared irregular | 0 | 19 | 30.1 | 23.1 | 19.9 | 20.9 / 1.7 |
| linear n=2 | few distinct | 0 | 3 | 27.7 | 22.8 | 19.6 | 17.9 / 1.7 |
| linear n=2 | fully irregular | 3800 | 3800 | 30.2 | 27.8 | 24.9 | 20.9 / 1.7 |
| linear n=6 | balanced | 0 | 1 | 77.0 | 83.1 | 68.8 | 31.7 / 2.7 |
| linear n=6 | shared irregular | 0 | 19 | 105.5 | 85.1 | 69.1 | 79.2 / 2.7 |
| linear n=6 | few distinct | 0 | 3 | 97.9 | 82.9 | 69.0 | 63.4 / 2.7 |
| linear n=6 | fully irregular | 3800 | 3800 | 105.4 | 99.4 | 84.1 | 79.2 / 2.7 |
| augmented CINT k=3 m=3 | balanced | 0 | 1 | 59.0 | 58.6 | 53.3 | 31.5 / 13.1 |
| augmented CINT k=3 m=3 | shared irregular | 0 | 19 | 88.3 | 59.0 | 53.1 | 79.0 / 13.1 |
| augmented CINT k=3 m=3 | fully irregular | 3800 | 3800 | 89.0 | 72.4 | 66.9 | 79.0 / 13.1 |
| augmented DRIFT k=2 m=4 | any | 3800 | 3800 | 99 to 103 | 85.7 to 85.9 | 81.9 to 82.9 | 93.5 / 29.0 |
| state-dependent n=2 | balanced | 7600 | 7600 | 49.5 | 45.5 | 40.8 | 40.5 / 9.1 |
| state-dependent n=2 | irregular | 8495 | 8495 | 54.1 | 50.0 | 45.3 | 44.9 / 10.1 |
| state-dependent n=6 | balanced | 7600 | 7600 | 192.5 | 166.6 | 167.1 | 191 / 18.7 |
| state-dependent n=6 | irregular | 8495 | 8495 | 212.5 | 186.4 | 183.3 | 213 / 20.8 |

GC time is 0% of every adjoint gradient in the second run (it was 6 to 40%).
The primal moved too, from the `Val` sites: linear n=2 balanced 7.0 to 5.5 ms,
1.4 to 0.5 MB. The 6-latent balanced regression in the first run was the
uncached packed LU and is gone in the second.

Every timing here is the minimum of seven repeats on a shared machine; treat
differences under about 5% as noise. The counts are not noisy.

## What remains, with measured shares

1. **Allocation left in the reverse pass.** After the scratch work the
   remaining sites are the state-dependent group replay
   (`_ctsem_complex_group_pullback!`, `_dual_state_derivative`, `_partial1`:
   tens of thousands of small allocations per gradient on augmented-drift
   models, 29 MB against 2.7 MB for a linear model of the same size) and the
   per-subject `:init` entry. GC is 0% of the adjoint in every measured cell
   now, so this matters for threads, not for serial time.
2. **`rev_update` at 13 to 29%** is now the largest single family on linear
   models. It has not been looked at.
3. **Lyapunov pullback for state-dependent models** (25% of the 6-latent
   state-dependent adjoint). With an orthonormal packed basis,
   `O(A') = O(A)'`, so the pullback could reuse the forward factorisation.
4. **Adaptive Padé degree in `my_exp!`** (always 13 today): the forward
   exponential is 4 to 6% of the adjoint wherever the cache misses.
5. **Augmented block exponential.** With `m` static carriers the Padé runs
   on `n = k + m`, but the carrier rows of `JAx` are zero, so
   `exp([J11 J12; 0 0] h) = [exp(J11 h)  φ1(J11 h) J12 h; 0  I]`, and the
   whole Padé recurrence (and the Fréchet one) can be carried in `k x k` and
   `k x m` blocks: `(k^3 + k^2 m)` per product against `(k + m)^3`, nine times
   less at k = 2, m = 4. It only pays where the cache misses, which for
   augmented models means individually varying drift or irregular intervals,
   and there `exp + frechet` is about 15% of the adjoint after this branch.
   Real but bounded; the structured Padé with its own Fréchet variant is
   roughly a day.

## Substepping: the mesh policy

The defect `r = f(x_1) - f(x_0) - J (x_1 - x_0)`, evaluated after a step from
`x_0` to `x_1` with the Jacobian `J` the step used, is the nonlinearity the
step actually encountered, per state, at the cost of one extra predict-group
evaluation and no exponential. It is identically zero for a linear model, so a
criterion built on it never substeps one, and it says which states are
nonlinear without any per-model configuration. Relative to the step's own
motion, `‖r‖ / max(‖x_1 - x_0‖, tol)` is a dimensionless indicator; a step is
refined by halving while the indicator exceeds a tolerance, up to a cap.

The mesh has to be smooth in θ for the optimiser and the Hessian, and it has
to be a stable shape for the tape. So the step count per row is **decided,
then frozen**: computed at the start of an optimisation pass from the
starting values, held fixed through the pass (line searches and gradients see
a fixed integrator), and recomputed between passes. ctsem already runs more
than one pass (prior warm-up, then likelihood), so the second pass gets a mesh
from a sensible point. A pass whose recomputed mesh differs from the one it
used reports how many rows changed, which is also the answer to "did
substepping do anything for this model": all zeros means no. `maxtimestep`
stays as a ceiling on the step, so the manual control still works and the
automatic mesh can only add steps below it.

Nothing in the reverse pass changes: the tape already records the substep
count each row actually took. The state-explicit path
(`intoverstates = FALSE`) is the natural oracle for the tolerance: it is the
model itself given the discretisation, one innovation per substep, and its
disagreement with the EKF fit on the same data is the linearisation error the
tolerance is meant to bound. Not implemented on this branch; the design is
here so that whoever implements it does not have to rediscover the smoothness
constraint.

## Not doing

- A singular-safe discretisation (Van Loan block) for exactly-zero drift.
  Decided: the `-1e-5` idiom stays.
- Higher-order exponential integrators or Magnus for the mean: the mean is
  already second order and the covariance recursion, not the mean, is the
  approximation.
- A spectral once-per-evaluation route for linear models on fully irregular
  intervals: real but fragile near defective drifts (critically damped
  oscillators); the kernel above takes the same case from 17% to 8%.
