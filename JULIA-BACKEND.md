# The Julia backend for ctsem

This file is the maintainer reference for the julia backend and is kept up to
date with the code. The vignettes are the user documentation. The report at
`../juliaBackendReport/julia-backend-changes.qmd` is frozen as a dated record
as of 2026-09-02 and is not updated; read this file or the vignettes instead.

`ctFit(..., backend='julia')` fits with a Julia extended-Kalman-filter engine
instead of Stan's generated C++, using a hand-written reverse-mode adjoint for
the gradient. This document records what it can do, where it deliberately
departs from the Stan backend, and how it is deployed.

A third, hand-written C++ engine was built and benchmarked alongside it; it
lives on the `cppBackend` branch, along with its own report (`CPP-BACKEND.md`)
and its Stan parity suite. Once the speed comparison levelled out, the Julia
engine was the more maintainable of the two, so this branch carries the Stan and
Julia backends only. Much of the architecture below was developed with both
engines in play, which is why it is written to be engine-agnostic — that shape
is worth keeping even with one engine, because the alternative is duplication
with Stan.

---

## Getting started

```r
ctJuliaInstall()
ctFit(data, model, backend = 'julia')
```

`ctJuliaInstall()` supplies whatever is missing — the `JuliaConnectoR` bridge
package, a Julia installation, and the vendored engine's Julia dependencies —
and skips whatever is not. Each step that installs something asks first, naming
the download and where it goes. Nothing needs a restart, and there is no
`JULIA_BINDIR` to set: a package installed into the running session's library is
visible immediately, and a Julia that ctsem installs lands under
`R_user_dir("ctsem", "data")`, which every later session searches.

`ctFit(backend='julia')` on a machine with none of this offers the same setup at
the point of failure rather than erroring, before it prepares any data — so the
first thing a user does with the backend is the thing they meant to do.

The three steps used to be four manual ones across two R sessions: install
`JuliaConnectoR`, install Julia, set `JULIA_BINDIR`, restart R, call
`ctJuliaSetup()`. Each was reported only as the error text of whichever call hit
it first.

Two things are worth knowing about the mechanism:

- **Nothing reaches the network without consent, which is what keeps this
  CRAN-legal and `R CMD check` offline.** An interactive session is asked; a
  non-interactive one *declines by default* and says how to consent in advance
  (`ctJuliaInstall(agree = TRUE)`, or `CTSEM_JULIA_AGREE=yes`). Refusal is the
  default rather than the fallback, so a scripted or automated run never
  downloads a quarter of a gigabyte by surprise.

- **The Julia version is pinned, not resolved at run time.** That is what lets
  the archive's sha256 ship inside the package and be checked before anything is
  unpacked, and the pin is the version the vendored `Manifest.toml` was resolved
  under, so the engine's dependencies instantiate as tested rather than being
  re-resolved on a user's machine. `ctJuliaInstall(version=)` overrides it, with
  a message saying the download went unverified.

Julia is downloaded only when none can be found. The search runs `JULIA_BINDIR`,
then a Julia ctsem installed itself, then `PATH`, then juliaup's own
installations — that last one because R started from a launcher rather than a
shell does not inherit juliaup's `PATH` entry, which is the normal case on
Windows. `ctJuliaStatus()` reports which one was chosen and what version it is,
and installs nothing, so it still works on the machine where nothing else does.

---

## The one thing that keeps the backends honest

Every backend consumes the *same canonical model specification*:
`.ctJuliaPrepare()` turns a `ctModel` into one row per model-matrix cell, with
each cell's transform rendered as an arithmetic expression string. The engine
interprets those strings at run time; it never re-derives what a parameter
means.

That is the load-bearing decision. The expensive, error-prone part of ctsem is
turning a `ctModel` into a canonical, augmented parameter table, and re-deriving
it per backend is exactly the seam where backends drift apart. Everything else
follows the same principle:

| shared by | mechanism |
| --- | --- |
| model preparation | `.ctJuliaPrepare()` — one parameter table |
| system-matrix summaries | `.ctSummaryMatricesFromArrays()`, factored out of `ctSummaryMatrices.ctStanFit()` |
| prediction output | `.ctKalmanArrayAssemble()`, factored out of `ctKalmanArray()` |
| uncertainty | `ctOptimComputeUncertainty()` — the Stan path's own, run at the end of `ctFit()` as `stanoptimis()` runs it |
| posterior predictive | `ctPostPredData()` / `ctFitCovCheck()` via `.ctFit*()` accessors |

Where a function needed backend-specific data it now reaches through an
accessor (`.ctFitModelObject`, `.ctFitLongData`, `.ctFitIdMap`,
`.ctFitObservedY`, `.ctFitRowSubject`, `.ctFitRowTime`,
`.ctFitObservedRowLoglik`, `.ctFitReplaceData`) rather than into `standata`.

---

## What works

- **Fitting**: `ctFit(backend='julia')`, continuous and discrete time,
  individually varying parameters (`intoverpop`), TI predictors, TD predictors,
  state-dependent (nonlinear) model matrices, `priors=TRUE`.
- **Threading**: the subject loop splits into contiguous chunks over
  `Threads.@spawn`, each owning its adjoint workspace — never indexed by
  `threadid()`, because a task can migrate between threads at any yield point.
  2.0–4.3x on six threads. `ctJuliaSetup(threads=n)` must run *before* the Julia
  session exists; `ctFit(cores=n)` caps it per fit.
- **Uncertainty**: `hessian`, `surrogate`, `is`, `opg`, `sandwich`,
  `bootstrap` — **run as part of fitting**, as `stanoptimis()` runs them for an
  optimized Stan fit, and controlled by the same `optimcontrol` names
  (`uncertainty`, `uncertaintyDraws`, `finishsamples`, `uncertaintyControl`,
  `estonly`). `ctOptimUncertainty()` re-runs it with different settings. The
  Hessian is **exact**, not finite-differenced — see below. Not
  `fullbootstrap`, which re-optimises each resample and so needs the model
  rebuilt rather than re-evaluated.
- **Cross validation**: `ctLOO()`, including `subjectwise`, `leaveOutN`,
  `keepfirstobs`, `refit=FALSE` and `casewiseApproximation`. Not
  `parallelFolds`, which is ignored: the engine threads its own subject loop,
  so each fold already uses every core.
- **Per-subject scores**: the adjoint already computes them; the rows sum to the
  full gradient to 1e-16.
- **Summaries**: `summary()`, `ctSummaryMatrices()`, `ctExtract()`,
  `ctDiscretePars()`, `ctSubjectPars()`.
- **Prediction**: `ctKalmanArray()`, `ctPredict()` / `ctKalman()`,
  `ctPredictTIP()`, `ctResiduals()`, `ctACFresiduals()`,
  `ctExtract(subjectMatrices=TRUE)`, `plot()`.
- **Generation**: `ctGenerateFromFit()`, and with it `ctPostPredData()`,
  `ctPostPredPlots()`, `ctFitCovCheck()`.

**Not supported**: `backend='julia'` refuses, rather than silently ignoring,
variational Bayes (`vb=TRUE`), data generation through the fit call
(`gendata=TRUE`), Stan compilation controls (`stanmodeltext`, `compileArgs`,
`forcerecompile`), manifest types beyond continuous/binary/ordinal/count/censored
(`manifesttype` outside `0:4`), and `intoverstates=FALSE` together with
`intoverpop='laplace'`. Also unsupported, without a hard refusal:
`fullbootstrap` uncertainty, `summary(priorcheck=)`
(accepted and ignored — it compares posteriors against the *Stan* model's
prior block), the prior/posterior-density, trace and interval panels of
`plot()`, and multi-start/restart robustness in the optimizer — both
`ctsem_optimize` and Stan's `stanoptimis` can walk to `|raw| ~ 1e4` on a weakly
identified model, and only Stan's has any defence against it. That last one
deserves more attention now that `ctLOO()` works, because every fold is an
independent re-optimisation: on a model where withholding a fold leaves a
parameter poorly determined, two folds — or two backends — can converge to
raw parameters several units apart for a likelihood difference of a percent.

---

## Transformed-parameter summaries

`summary()`, `ctSummaryMatrices()` and `ctDiscretePars()` all read the same
thing: samples of the *model matrices* implied by a raw parameter vector, as
`pop_DRIFT`, `pop_DIFFUSIONcov` and so on, each an `[iteration, row, column]`
array. So the architecture is one engine primitive and a thin stack on it:

| layer | what it does |
| --- | --- |
| `ctsem_parameter_matrices` | materialize every model matrix for an `npar x nsamples` block of raw vectors |
| `ctBackendParMatrices()` | the same, named and reshaped, for one vector |
| `ctExtract()` | the same over the posterior, as `pop_*` arrays |
| `ctSummaryMatrices()` | collapse those arrays — `.ctSummaryMatricesFromArrays()`, shared with Stan |
| `summary()` | fixed effects plus the system matrices |

The primitive runs **inside the engine**. The engine already materializes every
matrix from the raw vector on its way to a log likelihood, so asking it for that
same materialization is the only way to guarantee a summary reports what the
likelihood used; a second, R-side implementation of ctsem's transforms is
exactly the duplication this design exists to avoid. The derived matrices
(`DIFFUSIONcov`, `MANIFESTcov`, `T0cov`, `asymDIFFUSIONcov`, `asymCINT`) are
computed there too.

The batch shape is not incidental. The Julia bridge marshals a numeric array in
one transfer but a *list* element by element — a 344x difference measured
here — so a 200-draw posterior is one call, not 200.

Three things are worth knowing:

- **Verification is against Stan's own constrain step**, not against a second
  fit: `stan_constrainsamples()` at a fixed raw vector, on a model with an
  `intoverpop` augmentation, TI predictors and a state-dependent DRIFT. That
  takes the optimizer out of the comparison entirely.

- **`pop_T0VAR` differs from Stan, and it is a parameterisation difference, not
  a disagreement.** Stan computes `T0cov = sdcovsqrt2cov(T0VAR)` and *then*
  rescales `T0cov`'s indvarying-T0MEANS rows and columns by the parameter's
  multiplier and meanscale, leaving `T0VAR` unscaled. The engine folds that
  scale into `T0VAR`, so its `T0VAR` is the one whose `sdcovsqrt2cov` actually
  equals the reported `T0cov`. Both give an identical `T0cov` — which is what
  summaries report, and `summary()` drops `T0VAR` from the system matrices table
  for exactly this reason.

- **State-dependent cells are reported as conditional, not as constants.** They
  are functions of the latent state, so no single number describes them. They
  are evaluated at a state — T0MEANS by default, `state=` for anything else —
  and `attr(ctBackendParMatrices(fit), 'stateDependent')` names them, as does a
  note in `summary()`. This applies to any model with random effects, not only
  to explicitly nonlinear ones, because ctsem implements an individually varying
  parameter as a state dependence on an augmented carrier state.

Intervals appear only when earned: a fit built with `optimcontrol$estonly=TRUE`
has one "sample", and the interval columns are omitted rather than filled with a
zero-width interval that would read as certainty. An ordinary fit has 1000
draws, because `ctFit()` finishes with `ctOptimUncertainty()`.

---

## What `summary()` reports, and why each section is where it is

The sections are the same as an optimized Stan fit's, in the same order, with
the same names — `residCovStd`, `rawpopcorr`, `tipreds`, `parmatrices`, `popsd`,
`popmeans`, `logposterior`, `loglik` — because a script that reads
`summary(fit)$tipreds` should not have to know which backend produced the fit.
Getting there needed four things beyond the `pop_*` arrays above, and three of
them are the same trick: **displace the raw vector and read the population cell
back through the engine**, so the transform is applied once, by the code that
owns it, rather than reimplemented in R.

- **Which cell is a parameter's population value.** For most parameters it is
  the first cell they occupy. An `intoverpop` random effect is the exception,
  and not a rare one: the parameter's own cell is the carrier state's
  `T0MEANS[k]`, which holds the *raw* value and none of the transform, while the
  transform lives in whichever matrix reads `state[k]`. Reading the carrier cell
  reported a random-effects `CINT` parameter ten times too small, its `10*param`
  transform missing. Stan resolves the same ambiguity the same way, following
  its `pr2` reference before applying `tform`.

- **`popsd`** is the standard deviation of the *transformed* parameter over the
  population, which for a nonlinear transform is not the transform of the
  standard deviation. Stan draws 5000 subjects from the raw population
  distribution, pushes each through the transform and takes the sd. This does
  the same integral by 5-node Gauss-Hermite quadrature over the same
  distribution: exact for a linear transform, and steadier than a random cloud
  for a nonlinear one. All the varying parameters are displaced together on a
  shared node, so it costs one engine call per node rather than per parameter —
  only each parameter's own marginal spread is read, so the correlation a shared
  node induces between them never enters an answer.

- **`tipreds`** reports the effect on the *transformed* parameter, not the raw
  coefficient, which is what Stan's `linearTIPREDEFFECT` is: the transform
  evaluated a hundredth of an effect either side of the raw mean, differenced,
  and scaled back up. One predictor at a time, so it is two engine calls per
  predictor rather than two per effect.

- **`rawpopcorr`** is the one that needs no displacement. It is a correlation
  between *raw* parameters, so the state scaling folded into the sd transform
  cancels and `cov2cor(T0cov)` is already the reported quantity.

`residCovStd` comes from the filter rather than the transforms: the prior errors
at the estimate, standardised by the observed covariance. `ctFit()` caches the
filter output on the fit (`fit$kalman`) exactly as the Stan path caches
`stanfit$kalman`, so summarising does not repeat a whole filter pass.

**What a summary costs.** Nothing, now: about half a second on the model
above, against Stan's 0.8. It used to be sixteen seconds, because every call
re-materialized every model matrix for every draw. Stan was faster for a
structural reason rather than a clever one — it constrains its draws *once*,
at fit time, and `summary()` reads the cached `transformedpars` — and the fix
was to do the same. `ctFit()` now stores the constrained draws on the fit, and
`ctOptimUncertainty()` refreshes them when it changes the posterior; the cache
carries the draws it was built from, so it is checked rather than trusted.

That moves the cost rather than removing it: the constrain step is about
sixteen seconds at fit time on this model, paid once. Seven of its eight engine
calls are displaced reads — the five quadrature nodes behind `popsd` and the
two linearisation steps behind `tipreds` — and each of those materializes
*every* matrix, including the covariance factorizations and the solves behind
`asymDIFFUSIONcov` and `asymCINT`, when all it reads back is a handful of
scalar cells. An entry point that returns only the requested cells is the
obvious next reduction; it has not been measured against, so the size of the
win is unknown.

Measured against a Stan fit of the same model — five random effects, a TI
predictor, `intoverpop` — `residCovStd` agrees to the printed digits, `popsd`
and `popmeans` to ~1%, `rawpopcorr` to ~0.02, and `ctSubjectPars()` to 1e-4. The
residual differences are the two optimizers' own, plus Stan's Monte Carlo where
this uses quadrature.

Two consequences of fitting with uncertainty are worth stating, because both
are visible to a user who is only comparing backends.

**It costs.** The Hessian is 2 gradient evaluations per parameter, and it now
runs on every fit: the tutorial's 28-parameter individual-differences model went
from ~30 s to ~106 s. That is the same trade the Stan backend has always made,
and `optimcontrol$estonly=TRUE` opts out of it.

**Everything downstream now runs on the posterior**, where it used to run on a
single point estimate — and a draw from a normal approximation can land
somewhere the filter cannot go, e.g. a prior covariance the smoother's solve
finds singular. `ctExtract()` therefore drops inadmissable draws and reports the
proportion, which is what `stan_constrainsamples()` has always done; without
that, one bad draw in two hundred took the whole call with it.

---

## The Hessian, and cross validation

Two things the Stan backend does that this one reached for last, and that turned
out to need opposite kinds of work: one was a gap in the *engine*, the other a
gap in the *translation*.

### The Hessian is exact

`ctOptimComputeUncertainty()` gets its Hessian by central-differencing the
gradient: `2 * npar` reverse sweeps, accurate to about the square root of
machine precision, and only if the step suits the parameter's scale "—" which
one global step cannot do for a vector mixing log standard deviations with
unconstrained correlations. The engine can do better, because it can
differentiate its own gradient: `ctsem_hessian` runs `ForwardDiff` in forward
mode *over the reverse-mode adjoint*.

Forward-over-reverse rather than forward-over-forward because the reverse pass
is where the engine's work already is. `ForwardDiff.hessian` of the primal costs
`O(npar^2 / chunksize)` passes; differentiating the adjoint costs
`O(npar / chunksize)`, each one a traced forward sweep plus its reverse.

| | finite difference | forward-over-reverse |
| --- | --- | --- |
| relative error vs `ForwardDiff.hessian` | 2e-6 | **2e-16** |
| sweeps | `2 * npar` | `ceil(npar / chunksize)` |
| 28-parameter model, warm | 0.78 s | **0.08 s** |
| 28-parameter model, first call | 0.78 s | 10.6 s |

That last row is the one to know about. The first call on a given model spends
about ten seconds in Julia compiling the whole reverse pass for dual numbers,
and it is paid once per model shape per session "—" so a single fit in a fresh
session is slower, and everything after it is an order of magnitude faster. It
is announced (`Computing exact Hessian`) rather than left as a mysterious pause,
and `ctOptimUncertainty(control = list(analyticHessian = FALSE))` returns to the
finite difference. If the engine cannot differentiate at a point, the fallback
is automatic and warned about, because a worse covariance beats no fit.

**Making the reverse pass differentiable took four changes**, all of the same
kind: things that were `Float64` because nothing had ever asked them not to be.
The tape's records held the observed data, the TD predictors and the time step
as `Float64` beside `T`-typed siblings; the group-replay scratch in the dual
context was `Float64` with a comment explaining that the adjoint is only ever
entered with doubles; and the Fréchet block called `Base.exp`, which has no
method for a matrix of duals, where the engine's own `_ctsem_expm` already
dispatched correctly. None of it changes a `Float64` result "—" the engine's
815 existing tests pass unchanged "—" and `test_hessian.jl` adds the identity
that matters, plus a check that an ordinary gradient still works afterwards
(the workspace cache is keyed on the scalar type, so a Hessian call swaps it out
and the next gradient swaps it back).

### Cross validation

`ctLOO()` needed nothing from the engine. The Stan path withholds rows by
zeroing `standata$dokalmanrows`, refits with `stanoptimis`, and reads `llrow`
out of `constrain_pars`; none of those exist here, but each has an exact
equivalent, and only the middle one needed a function extracted:

| Stan | julia |
| --- | --- |
| `dokalmanrows[i] <- 0` | set row `i`'s manifests to `NA` and re-prepare |
| `stanoptimis(estonly=TRUE)` | `.ctJuliaOptimise()`, the call `ctFit()` makes |
| `constrain_pars(...)$llrow` | the filter's own `llrow`, a byproduct of the forward pass |

Withholding at the *data* level is what makes this short, and it is also
stricter than a flag: a withheld row cannot leak into the likelihood, because
the engine never receives its value. It is the same mechanism `removeObs`
already uses for prediction, so it is a tested path rather than a new one.

The identity worth asserting is the one with no optimizer in it: with
`refit=FALSE` every fold scores the same parameters against the same full data,
so assembling the out-of-sample vector fold by fold must reconstruct the
in-sample one exactly. Any error in which rows a fold owns shows up there with
nothing to hide behind. The comparison against Stan is made at `refit=FALSE`
for the same reason "—" with refitting, each fold is an independent
optimisation, and on a fold that withholds a third of the data the two
optimizers land up to 6 raw units apart in one parameter for a ~2% likelihood
difference. That is the known multi-start weakness, not a disagreement about
cross validation.

---

## Prediction: ctKalman, ctPredict, and subject parameters

Those functions read four arrays — prior, filtered and smoothed states and
observations for every data row — and every one is a byproduct of the forward
pass the likelihood already makes. So prediction rides on that pass rather than
adding a second filter "for prediction", which is precisely how the prediction
output and the likelihood would drift apart.

The mechanism is the recorder hook the filter already carries for the adjoint
tape. A trace is whatever implements the `_record_*!` hooks: the adjoint tape
implements the ones it needs, `CTSEMKalmanTrace` implements the ones it needs,
and every hook the other does not want resolves to an inlined no-op. Dispatch
rather than branching means the primal and adjoint paths compile to the code
they did before any of this existed.

Three hooks are new, because the adjoint has no use for what they capture. They
sit in the filter *loop*, not inside the measurement update, so they see the row
index and fire on a fully missing row — which the update returns early from, and
which is exactly where a recorder one level down would silently skip.

Smoothing is a genuinely separate backward RTS pass: it needs a subject's last
row before it can produce the first, so it is written as one rather than
disguised as part of the forward loop.

**Subject-level parameters** fall out of the same pass. ctsem represents an
individually varying parameter as an augmented latent state with no drift and no
diffusion, so a subject's value for it *is* its smoothed t0 estimate; a
subject's matrices are the parameter vector its filter ended with, T0MEANS
replaced by that state. This is Stan's construction, comment included ("t0means
updated, other pars as per final time point").

`ctSubjectPars()` reads those matrices, and reads each parameter from the same
population cell the fixed-effects summary uses — so a random-effects `CINT`
parameter comes back from `subj_CINT` with its transform applied, not from its
raw carrier state. It agrees with Stan's to 1e-4 on a five-random-effect model.

`removeObs` withholds observations from the filter but not from the report,
which is the point of it: what comes back is a prediction next to the
observations it was not given. `ctResiduals()` and `ctACFresiduals()` come along
for free, since they are built on `ctKalmanArray(standardisederrors=TRUE)`;
their standardised residuals have sd 1.00 at the estimate of a correctly
specified model, which is an independent check on the whole chain that no
comparison against Stan provides.

`ctPredictTIP()` predicts at chosen covariate values by building a dataset of
pseudo-subjects, one per value, and asking the fitted model for its expectation
on it. That is a *data* operation, so once a fit can be re-prepared against a new
data frame the whole function follows; its dynamics panels route through
`ctDiscretePars()` with per-subject matrices, so they exercise the
subject-parameter path as well.

---

## Three places this improves on Stan rather than reproducing it

None of these touches the likelihood or the gradient — all three change only
what is *reported* — and each is verified by an identity rather than by a
comparison against the code that produced the numbers.

**1. The measurement model is re-evaluated at the updated state** before the
filtered observation estimate is recorded. Stan applies LAMBDA and MANIFESTMEANS
as the *prior* state left them, and its own comment on that block reads "these
could be improved by recomputing all state dependent pars, error covariances
etc. at each step". This is the ordinary case rather than an exotic one: ctsem
represents an individually varying parameter as an augmented latent state, and
MANIFESTMEANS is individually varying *by default*, so the measurement intercept
is a latent state whose pre-update value Stan reports. On a two-process model
with ctsem's defaults the reported `yupd`/`ysmooth` move by ~1 unit.

The test: when the measurement equation is linear in the augmented state —
which every intoverpop model's is — `y = Jy x` exactly, at prior, filtered *and*
smoothed. That holds now and did not before.

**2. The interval transition composes the substeps.** The filter propagates
`x <- A_s x + b_s` at each bounded substep, so the interval Jacobian is
`A_S ... A_1`. Stan instead recomputes `exp(JAx dt)` from the last substep's
Jacobian. For a state-independent JAx the two agree exactly; for a
state-dependent one only the product is the derivative of what was actually
computed. Composing is also cheaper — one matrix multiply per substep rather
than an extra exponential.

**3. The transition includes the TD impulse Jacobian.** An impulse sits between
the previous row's posterior and this row's prior and maps the covariance
through `Jtd`, so it belongs in the transition the smoother uses. Stan saves
only the exponential. `Jtd` is the identity in the common case; the
non-identity case, which a state-dependent TDPREDEFFECT produces, is tested in
the engine suite.

**Consequence:** `ctKalman()`/`ctPredict()` output for a model with individually
varying measurement parameters will *not* match Stan's, by design. A test
asserts exactly which parts diverge and which — latent states, likelihood, prior
estimates — still agree.

One further deliberate difference, in the filter itself: for the **first row of
each subject**, the manifest covariance is built from MANIFESTVAR *after* the
update-group transforms rather than before, so a state-dependent MANIFESTVAR is
not read before its transform has ever written it.

The first row also runs all three transform groups -- predict, td, update -- in
the order the main loop runs them, even though it has no prediction interval. A
group supplies values as well as a prediction: PARS is in the predict group, and
an update-group or td-group cell may be written by a transform that reads a PARS
cell. Their contexts carry a zero interval, which no transform can observe --
`generate_complex_transform_string` substitutes only states, PARS and the model
matrices, so no transform expression can name the interval or the time.

---

## Data generation

Generation is the one thing here that is not a passive recorder: it *changes*
what the filter consumes. Each row's observation is drawn from its own prior
predictive as the filter reaches it, and the filter then carries on as though
that draw had been read from the data — so the state it propagates is
conditioned on the drawn history, not the real one. That is exactly what makes
the result a draw from the model rather than a sequence of independent one-step
predictions, and it is also the easy thing to get wrong: a simulator that
predicted each row from the *real* history would produce data that looks
perfectly plausible and is a draw from nothing.

The draw is taken from the innovation covariance the update is about to use,
*after* it has already been factorized, so the drawn innovation is `L z` by
construction and cannot come apart from the covariance the filter then
conditions on.

The standard normals are drawn in R rather than by an engine RNG, so
`set.seed()` means what a user expects.

Two properties are asserted rather than assumed:

- **The defining identity.** Re-running the ordinary likelihood on the generated
  dataset reproduces the likelihood reported while generating it. Nothing that
  drew from the wrong covariance, or conditioned on the wrong history, satisfies
  that.
- **Calibration.** Over 200 generated datasets the observed data's log
  likelihood sits at the 54th percentile — an unremarkable draw from the fitted
  model, which is what it should be at the maximum.

Missingness is preserved: an entry that was not observed is returned NA rather
than invented, because a posterior predictive check compares against the
observations that exist.

---

## Discrete time

`type='dt'` models work, and did not need a second filter.

A discrete model's DRIFT, CINT and DIFFUSION *are* the one-step quantities, so
every expensive piece of the continuous form collapses: the transition is `JAx`,
the process noise is the diffusion covariance, and the intercept is the local
affine offset with no solve around it. Both the forward and the reverse pass
come out **shorter** than their continuous counterparts, and everything
downstream — the measurement update, the smoother, subject parameters,
prediction, generation — is untouched. The branches are separate functions
rather than a test inside the recursion, which would put one in every step of
it.

The asymptotic forms follow: `(I - A) x = c` for the intercept, and the discrete
Lyapunov equation `X = A X A' + Q` for the asymptotic diffusion. `dtDRIFT` is
absent from the summary, because there is nothing to discretise.

The property that separates discrete from continuous is asserted directly:
stretching the recorded times changes nothing, because each row advances exactly
one step, while a continuous reading of the same data does change.

---

## How it is verified

| suite | what it checks |
| --- | --- |
| `tests/testthat/test-stan-julia-parity.R` | likelihood and gradient against Stan across model shapes |
| `test-backend-summary.R` | `pop_*` arrays against `stan_constrainsamples()` at a fixed raw vector; that a default fit carries uncertainty and reports the Stan sections |
| `test-backend-kalman.R` | per-row prior/filtered/smoothed output, subject matrices, `ctPredict`, `ctPredictTIP` |
| `test-backend-generate.R` | the generation identity and calibration |
| `test-backend-discretetime.R` | discrete time against Stan, and the interval-independence identity |
| `test-backend-uncertainty.R` | Hessian standard errors against Stan's; the exact Hessian against the finite difference it replaces |
| `test-backend-loo.R` | `ctLOO()` folds, the no-refit identity, and agreement with Stan |
| `test-backend-priors-scores.R` | `priors=TRUE` against Stan; scores sum to the gradient |
| `test-julia-engine-vendored.R` | `ctJuliaSetup()` from the vendored copy, offline |
| `test-julia-install.R` | download URLs, archive unpacking, and that consent is refused rather than assumed when there is nobody to ask |
| `inst/julia/ContinuousTimeSEM/test/` | the engine's own suite, including the adjoint against ForwardDiff and `test_hessian.jl` |

Two habits worth keeping. Comparisons against Stan are made **at a fixed raw
parameter vector** rather than between two fits, so the optimizer is not part of
the comparison. And where this backend deliberately differs from Stan, the test
asserts the *identity the better version satisfies*, not merely that the numbers
changed.

The installed `ctsem` in a user's library will usually lag this working tree, so
install it to a throwaway library before running any of these:

```bash
R CMD INSTALL --library=<scratch> --no-multiarch --no-docs --no-byte-compile ctsem
```

and point `CTSEM_JULIA_PROJECT` at a Julia checkout if you want to test engine
changes before vendoring them.

---

## Deployment

The engine **is part of ctsem**: `inst/julia/ContinuousTimeSEM/` is its source,
edited in place like any other file in this repository. Originally it was
installed with `Pkg.add` from a private GitLab URL, which meant:

- **`backend='julia'` could not be installed by anyone outside one project.**
  An unauthenticated `git ls-remote` returned *HTTP Basic: Access denied*. It
  worked locally only because of an access token beside the checkout — a wall,
  not friction, and invisible from inside the group, which is why nothing caught
  it.
- **`revision` was a branch name**, so the "lock file" locked nothing.
- **The engine lived only in a second repository** that had to stay in lock-step
  with `R/ctJuliaBackend.R`'s parameter-table contract by hand.

The first fix copied the engine in and recorded the upstream commit in
`inst/julia/engine.json`, refreshed by a `tools/sync-julia-engine.sh` script.
That removed the credentials problem but kept the two-repository one, and added
a worse failure of its own: the cached engine project is keyed on that recorded
revision and only populated when empty, so editing the engine without also
running the sync script left everyone who had already used the backend running
new R code against their old cached engine, silently. An identifier that has to
be maintained by hand is an identifier that will eventually be wrong.

Both are now gone. There is one repository, no lock file, and no sync step:

| | originally | vendored + lock | now |
|---|---|---|---|
| install source | `Pkg.add` from a private URL | copy + `engine.json` | part of ctsem |
| credentials / network at setup | required | none | none |
| engine identity | a branch name | a recorded commit | a hash of the engine source |
| keeping the copy current | n/a | `tools/sync-julia-engine.sh`, by hand | nothing to do |
| dependency closure | 111 packages | **67** | **67** |
| fresh depot | 268 MB | **124 MB** | **124 MB** |
| user steps from nothing to a fit | 4, across 2 R sessions | **1** | **1, `ctJuliaInstall()`** |

Most of the dependency weight was DataFrames, used only as a row container in
the R interface, plus four dependencies with no call sites at all that were
loaded because the module opened with `@reexport using` over all of them.
`ekf_from_columns` replaced `ekf_from_data_frame`, taking the parameter table as
plain column vectors with sentinels (`0`, `NaN`, `""`) instead of `missing`;
DataFrames stays as a *test-only* dependency, because a `DataFrame` literal is
still the clearest way to write a small table in a test.

`ctJuliaSetup()` copies the engine into a writable project under
`R_user_dir("ctsem", "cache")` and instantiates it there — activating it in
place would fail on a read-only R library. If the manifest is unsatisfiable on
the user's Julia version it falls back to a fresh resolve rather than refusing
to run.

**That project directory is keyed on `.ctJuliaEngineVersion()`, a hash of the
engine's own source.** This is the piece that makes editing the engine safe:
change any file under `inst/julia/ContinuousTimeSEM/` and the key changes, so
the next `ctJuliaSetup()` builds a fresh project instead of reusing a stale one.
Nothing has to be remembered, which is the only property that survives contact
with actual use. `ctJuliaStatus()$engine` reports it.

### Working on the engine

Edit `inst/julia/ContinuousTimeSEM/` and commit it with the rest of ctsem.
Nothing else is required — no sync, no version bump.

```bash
# run the engine's own Julia suite against the tree you are editing
julia --project=inst/julia/ContinuousTimeSEM -e 'using Pkg; Pkg.test()'
```

`ctJuliaSetup(project = "<path>")` points ctsem at a different checkout, which
is useful for comparing against another copy but is not needed for ordinary
work.

### Two bridge behaviours to know about

- **JuliaConnectoR hangs — does not error — marshalling a zero-length vector**,
  in either direction. Optional parameter-table columns are Julia keyword
  arguments passed only when they have entries, and
  `ctsem_parameter_layout` reports a *count* of state-dependent cells with the
  cells themselves fetched separately, because a linear model has none.
- **A length-one R vector marshals as a scalar.** `.ctJuliaVector()` sends those
  as a list so they arrive as an `AbstractVector`; a single-subject objective —
  which `ctPredict()` routinely builds — otherwise hands the constructor an
  `Int` where it wants a vector.

### Spinning the engine out as a standalone Julia package

`inst/julia/ContinuousTimeSEM/` is laid out exactly as a **root-level Julia
package** — `Project.toml`, `Manifest.toml`, `src/`, `test/`, and nothing else.
That is deliberate, and it is what keeps the single-repository arrangement from
being a one-way door. Because the directory is a package root and its contents
live in ctsem's own git history, `git subtree` can publish it to a separate
remote whenever that becomes useful, and pull changes back:

```bash
git subtree split --prefix=inst/julia/ContinuousTimeSEM -b julia-engine
git subtree push --prefix=inst/julia/ContinuousTimeSEM <remote> main
git subtree pull --prefix=inst/julia/ContinuousTimeSEM <remote> main --squash
```

Someone can then fork or `Pkg.add` the engine on its own, with no R involved,
while ctsem development stays a single-repository affair. Nobody has to run any
of this to work on ctsem — it is available if the engine ever wants a life of
its own.

A **submodule** would be the wrong tool here, and it is worth saying why since
it is the more obvious reach: submodule contents are not part of the parent
repository's tree and are not included in an R package tarball, so a released
ctsem would ship an empty directory and users would be back to needing network
access and credentials — the exact problem this arrangement exists to avoid.
Subtree keeps the content in the tree while leaving the two-way link available.
