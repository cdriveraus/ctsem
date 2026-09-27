# dev/optimbench: the optimiser bench

The instrument every change to the julia optimisation flow is measured against
(`review/OPTIM-consolidation-plan-2026-09-25.md`, P0). A fixed set of problems,
fixed data, fixed starts, deterministic counts, a runner for dev1, and a
summary that says which machine produced every number. Build-ignored (`dev/`).

## The cell contract

A cell is one line of a grid file and runs as one R process:

    id  model  data  route  start  variant

plus the build it runs against (a label, see below).

- **model** names a problem in `cells.R`: a simulator and a model builder.
  Families: the stochopt regimes (`panel`, `panel5k`, `long`, `ordinal`,
  `nonlin`, `small`, `bigp`), the carefulfit/gaptol design (`cf_gaussian`,
  `cf_binary`, `cf_ordinal`, `cf_mixed`), the gated-gaps configs (`gA1`-`gA16`,
  `gB1`-`gB8`, `gC1`-`gC8`, `gD1`-`gD4`, `gN1`-`gN4`), AnomAuth (`anomS1`,
  `anomS2`), and three test fixtures (`acnonlin`, `mvmix`, `jflat`). Each says
  where it came from.
- **data** is the simulator's seed, or `cfg` for a problem that fixes its own.
  No simulator calls `ctGenerate` (its draw stream moves under unrelated
  commits). Data are simulated once into the store (`~/dev/ctsem-bench-data`)
  and read from there by every later cell and every build, and each result
  records the stored file's md5, so two builds can be shown to have seen the
  same rows.
- **route** is `intoverpop`: `augmented`, `laplace` or `auto`.
- **start** is `default:<seed>` -- `inits = NULL` with `set.seed(seed)`
  immediately before `ctFit`, the path a user takes -- or `zeros`, or
  `stored:<name>` for a vector in `starts.R` (the AnomAuth spurious maxima, the
  flat-transform start, each gated-gaps config's best-known point). Supplying
  `inits` switches the prior warm-up off, so a stored start is never the
  default path.
- **variant** names an entry of `variants.R`: `optimcontrol` and argument
  overrides, the factor under test. A variant that changes the objective (a
  prior, a floor) must say so with an `objective` tag, or its losses are
  measured against a different function.

Every fit uses `cores = 1`: at `cores > 1` the chunk tuner times candidates and
the estimate moves by ~1e-7 between runs.

### What a cell records

In `<results>/<id>.rds`, a list whose `row` is one row of scalars (the CSV
columns) and the rest the full record:

- the build (label, sha, ctsem version), the harness sha, machine, pid and
  Julia pid, and the 1-minute load at process start, fit start and fit end;
- the **contamination control**: seconds per value-and-gradient evaluation of
  this cell's own objective at raw zeros, timed inside Julia in three batches
  of about a quarter second each, after an untimed call that pays the
  compilation. No optimiser change can move it;
- the wall seconds of a warm-up fit, which is the timed fit itself run once
  first (same data, arguments, start and seed), so that no stage the timed fit
  reaches compiles inside it, and whether the two fits agreed (`warm_dx`,
  `warm_dll`); then of the fit, and of its optimiser, certification,
  uncertainty and Laplace-correction stages. A warm-up on a slice of the
  subjects capped at five iterations, what this was until 2026-09-26, left
  every later stage to compile in the timed fit (gA1: 6.9 s at one build,
  113 s at the next with identical counts);
- **per stage**, from thin wrappers the harness puts in the namespace (they
  forward every argument untouched): each `.ctJuliaOptimise` call with its
  kind (`warmup`, `main`, `resume`, `restart`), caller, iterations, objective
  and gradient calls, batch sizes and Newton steps; each engine optimisation
  run, including the stall-escape stages whose counts the stage's own result
  does not sum; each Hessian computed, and whether a reused one was only
  handed on. `warmup_ran` comes from this record; `warmup_claimed` is what
  `fit$optim$carefulfit` says, and the summary shows the two side by side;
- from the fit: `$optim` (without the trace, which is kept separately),
  certification status and gap, log likelihood and log posterior, the raw
  estimate and standard errors, the Laplace correction record and
  conditioning, the identifiability count, the engine's operation counts
  (`ctsem_opcounts()`), every warning and the first 200 messages;
- the end point **re-scored on a freshly built objective** (cold inner modes);
  for Laplace cells the 5-node quadrature there (`ctLaplaceCheck`), and the
  **exact reference**: per unit, `bench_probe_reference` (softcut 3.5, as the
  gaps note settled) when units have dimension at most 3, else importance
  sampling with two proposals (the nested configs, d = 11); plus the prior
  term, so the value is the penalised exact log likelihood the objective
  approximates. Each computed reference is checked against importance
  sampling on the three units where it departs most from Laplace. Seeds that
  reach the same point (raw within 1e-4) share one reference computation;
- `status`: `ok`, `evalonly`, `control`, `error` or `timeout`. The fit is
  capped at `BENCH_TIMEOUT` (7200 s) and the references at `BENCH_REFCAP`
  (5400 s); a cell killed from outside shows in `status.tsv` and the summary
  as `killed(exit 124)`.

### How cells are scored

`summarise.R` scores a cell by its penalised exact log likelihood where it has
one, else its re-scored objective, and compares it only with cells of the same
model, data, route, objective tag and score kind -- across every results
directory it is given, and for exact scores against the best-known value in
`references.csv` too. `loss` is the score minus that best (0 is the best
known), `dx` the largest raw distance from the best point.

`references.csv` holds, for each included gated-gaps and AnomAuth config, the
best penalised exact value any floor reached from any start in the gaps job's
sweeps, or any bench cell since where one beat it, and says which fit; the
stored point in `starts.R` moves with it. When a baseline's best exact score
for a config beats the file (by more than 1e-3, or 0.05 for an importance-
sampling reference), move both. The grid re-evaluates each at its stored point
(`variant = evalonly`), which checks the data, the model, the prior term and
the reference code together: at juliaFit `e2abf637` config A1 reproduced
-366.9738 to four decimals.

## Running it

### Once: the standing tree on dev1

    bash dev/optimbench/ship.sh --init

makes `~/dev/ctsem-bench`, a git clone of juliaFit from a bundle (dev1 cannot
reach this machine), with `src/` built for Linux: it generates the Stan
sources with dev1's rstantools and reuses the objects of an existing dev1 tree
whose generated `.h` files are byte-identical and whose `inst/stan` and
`inst/include` digest the same (the `.cc` are thin wrappers and prove
nothing). It installs to `~/dev/ctsemlib-bench-standing`.

### Each build: ship it

From any worktree of the ctsem repository:

    bash dev/optimbench/ship.sh <commit-or-branch> <label>

bundles only what dev1's clone does not have, fetches it there, checks the
commit out at `~/dev/ctsem-bench-wt/<label>` with the standing clone's `src/`
(no Stan rebuild unless the build changed `inst/stan` or `inst/include`), and
installs it to `~/dev/ctsemlib-bench-<label>` with a `BENCH_BUILD` file. A
label is bound to one sha; use `<branch>-<shortsha>`. The harness itself is
shipped separately and runs against any build:

    bash dev/optimbench/ship.sh --harness            # the optimbench branch tip

### Run a grid

    ssh cd-dev1 'bash ~/dev/ctsem-bench-harness/dev/optimbench/run_grid.sh \
      ~/dev/ctsem-bench-harness/dev/optimbench/grids/baseline.txt <label> 16'

returns at once; the runner detaches, writes its pid to
`~/dev/ctsem-bench-results/<label>/baseline/GRID.pid` and its log beside it,
first stores every cell's data and loads the build's engine once, then runs
the cells 16 at a time, each as its own process with one Julia thread. It runs
from a copy of the harness taken at launch, so shipping a new harness does not
touch a running grid. `GRID_DONE` appears when every cell has finished.
Launching the same grid again resumes it (finished cells are skipped). Stop it
with `stop_grid.sh <results dir>`, never `pkill`.

`grids/smoke.txt` checks a build and the harness end to end in minutes. To
make a grid for your own factor, copy one, change the `variant` column (adding
the variant to `variants.R`), and ship the harness.

### Summarise

    ssh cd-dev1 'cd ~/dev/ctsem-bench-harness/dev/optimbench && Rscript summarise.R \
      --out ~/dev/ctsem-bench-results/<label>/baseline/summary \
      ~/dev/ctsem-bench-results/<baseline label>/baseline ~/dev/ctsem-bench-results/<label>/baseline'

writes `summary.md` and `summary.csv`. The first directory is the baseline;
cells present in both are paired. The report prints the contamination check
first, then whether each fit's warm-up flag says what ran, then medians by
model, route and variant, the reference checks, and every cell.

Summaries that are kept go in the CT-SEM root repository under `review/bench/`
(the raw `.rds` gitignored there).

### Locally

A cell runs on Windows against a source tree, for correctness only:

    BENCH_TREE=<worktree> BENCH_DATA=<scratch>/data BENCH_LABEL=local \
      Rscript dev/optimbench/harness.R t1 cf_gaussian 1 laplace default:1 default <scratch>/out

Never time anything locally, and never edit `harness.R` while a cell runs from
it: Rscript reads a script as it executes it.

## Rules the numbers depend on

- Timing only from dev1, only when nothing else runs there, and every timing
  names its machine. Prefer the counts: iterations, objective and gradient
  calls, Hessians, operation counts.
- Read the contamination check before any timing. If the control moves by
  more than the effect, there is no timing result.
- Paired comparisons only: same cell id, same data md5, same start.
- A variant that changes the objective carries an `objective` tag.
- `fit$optim` flags are what the fit says; the stage record is what ran. When
  they disagree, the stage record is right and the flag is a defect.
