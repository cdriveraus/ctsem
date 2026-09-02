# Stan deprecation: where the seams are

The stan backend stays for now. It is what users know, and every existing
script, tutorial and the training material assume it. It will go once the julia
backend has carried users through a release or two. This file records, as they
are met, the places where the stan path costs something today and what removing
it would free, so that when the time comes the work is a list and not an
archaeology dig.

Nothing here is a recommendation to remove anything now. Reviews add rows;
nobody acts on them until Charles calls the deprecation.

## Rules for adding a row

- A row is a place in the code, not an opinion. Give `path:line` at the commit
  you read.
- Say what the stan-only code serves: a feature, a compatibility path, a
  data format, a class name users pattern-match on.
- Say what its removal would free: a dependency, a compile step, a duplicated
  function, a branch in a shared function, a confusing argument.
- Say what a user would notice. "Nothing" is a valid and useful answer.
- If the julia path is missing something the stan path has, that is a
  *blocker* row, not a removal row; mark it so. Deprecation cannot start while
  any blocker stands.

## Removal points

| location | serves | freed by removal | user notices | found by |
|---|---|---|---|---|
| `DESCRIPTION` LinkingTo rstan, StanHeaders, BH, RcppEigen, RcppParallel; Imports rstan, rstantools | compiling `inst/stan/*.stan` at install | the 20-minute compile, the mingw object-size limit that blocks `load_all()` with compilation, most of `src/` | install time only | survey |
| `configure`, `configure.win`, `R/stanmodels.R`, `src/` | rstantools generated model objects | the generated-code step and its Makevars | nothing | survey |
| `R/ctTIpredAuto.R` whole file (`scorecalc`, `ctTIauto`, `whichsubjectpars`) | stan-path per-subject scores, by re-fitting a stan model per subject over a PSOCK cluster; callers `R/stanoptimis.R:1020`, `R/ctLOO.R:80`, `R/ctOptimUncertainty.R:6,308` | the file, once the stan `scorecalc()`-based uncertainty, LOO and TI paths retire; the julia engine already returns per-subject scores directly (`ctsem_subject_gradients`, used by `R/ctBackendUncertainty.R:442-452` and `ctIdentify()`) | nothing, if the julia score path covers the same callers | J1.6b |
| seven stan log-prob-plus-gradient wrappers: `R/stanoptimis.R:407-417,502-512,596-611,630-643,1431-1487`, `R/ctOptimUncertainty.R:1093-1134` | the stan objective under optimisation and uncertainty | about 150 lines, and with them the single-core versus multi-core inconsistency where only the multi-core wrapper clips an infinite gradient and traps a NaN one | nothing | J1.3 F1, F3 |
| parsteps auto-freeing machinery, `R/stanoptimis.R:203-560` | stan-path staged freeing of parameters; untested, and its roxygen points at a `ctFitAuto` that does not exist | about 350 lines | only if julia has no equivalent staged fit; check whether `carefulfit` covers the need first | J1.3 F4 |
| `ctFit(vb=TRUE)` branch, `R/ctFit.R:1232-1246,1280` | Stan variational Bayes. It crashes on the default `fit=TRUE` because the branch never sets `rawest`/`rawposterior`, and it is silently ignored when `optimize=TRUE` | the branch and the `vb` argument | an error for a call that already errors | J1.4 F1 |
| `ctCheckFit()`, `ctFitMelt()`, `R/ctFitCovCheck.R:558-998` | stan-only generated-versus-observed diagnostic, about 440 lines including large commented-out former features | the code; `ctFitCovCheck()` and `ctPostPredPlots()` cover the purpose on both backends | `ctCheckFit` is exported, so it needs a deprecation message pointing at `ctFitCovCheck` | J1.5 F4, F13 |
| `ctChisqTest()`, `R/ctCompare.R:100-109`; `ctModelCoverage_check()`, `R/ctCoverageCheck.R` | stan-only, no backend check, no functional test | two exports | removal needs a deprecation cycle | J1.5 F7, F8 |
| `R/stanWplot.R`, `R/stan_unconstrainsamples.R`, and the `stan_reinitsf` call sites | Stan HMC plotting and re-initialisation helpers | the files | nothing for julia users | J1.4 |
| `R/ctStanParMatrices.R` | stan-only `expm(DRIFT*t)`, reachable only through `ctTIpredEffects()` | the file; `ctDiscreteParsDrift()` in `R/ctDiscretePars.R:317` is the shared version | nothing | J1.2 F8 |
| `R/ctData.R:88-149,197-206` and `R/ctFit.R:1046` | Stan's data convention: a missing TI predictor becomes 99999 for the generated model to impute, a missing TD predictor becomes 0 | the placeholder convention and the pre-fill. The julia path is handed the same pre-filled data and cannot read the sentinel, so it currently fits a covariate value of 99999 | nothing once julia handles missing values itself | J1.1 F1, V2 |
| `ctStanModelWriter` and its five nested closures, about 1,235 of `R/ctModelWriter.R`'s 2,008 lines | generating the Stan model text | those 1,235 lines. The rest of that file, despite its `ctStan*` names, is shared infrastructure the julia path calls directly, so the file splits cleanly rather than disappearing | nothing | J1.6a Q1 |
| `vignettes/hierarchicalmanual.rnw`, `vignettes/optim-uncertainty.qmd`, the README's Rtools section | documents written when stan was the only path | the stan-specific prose. The Rtools requirement goes with the compile step | documentation reads as julia-first | J4.1 F7, F8, F9 |
| `fit$stanfit$transformedparsfull` versus `fit$transformedpars` (`R/ctFit.R:1286`, `R/ctJuliaBackend.R:2472`) | two storage locations for constrained draws, one per backend | one location, or one accessor users are told to use instead of either | users who read the stan slot directly need a pointer; `ctExtract()` is that pointer today | J4.1 F5 |

## Blockers

Things the julia path does not yet do that the stan path does. Each must be
resolved or explicitly dropped before deprecation starts.

| capability | stan location | status | found by |
|---|---|---|---|
| plotting a sampled fit's draws | `R/ctStanPlotPost.R`, `R/plot.ctStanFit.R` | no function reads a sampled julia fit's `estimate$rawposterior`; `ctTracePlot` covers the optimiser trace only. Also undocumented: stan takes the median of draws as the point estimate and julia the mean | J1.4 F6, F9 |
| `fullbootstrap` uncertainty | `R/ctOptimUncertainty.R` | not supported on julia, because each resample needs the model rebuilt rather than re-evaluated | survey |
| `summary(priorcheck=)` | `R/summary.ctStanFit.R` | accepted and ignored on julia; it compares posteriors against the *Stan* model's prior block | survey |
| `ctEmpiricalBayesFit()` | `R/ctEmpiricalBayesFit.R:82,90-99` reads `fit$stanfit$rawest` and `fit$standata` | stan-only, and until this review said so nowhere. Needs a julia port or an explicit decision to retire it | J1.3 F5 |
| `ctTIpredEffects()` | `R/ctStanTIpredeffects.R:67-71` reads `fit$ctstanmodel` | fails on julia fits. `ctPredictTIP()` and `summary()` already report TI effects there, so this may want a pointer rather than a port | J1.2 F1 |
| analytic Jacobian of the drift for nonlinear models | `R/ctJacobian.R`, consumed by `R/ctJuliaBackend.R:615-618` to build the engine spec, and read back by every nonlinear helper through `ctBackendParMatrices()` | not a removal point at all: shared by both backends today. The engine has no Jacobian source of its own, so this R file stays until it does | J1.6b F4 |

## Resolved or reframed

| item | outcome | found by |
|---|---|---|
| HMC diagnostics parity | Reframed. The stan path never surfaced divergences, tree depth or step size at all: `get_sampler_params` is imported and never called. The julia path returns them. Julia is ahead here, and what remains is the plotting row above | J1.4 F5 |
| variational Bayes | Dropped as a blocker. The stan implementation is broken, so there is no parity to reach | J1.4 F1 |

## What to protect rather than delete

The julia engine's numerical conventions are pinned to Stan's generated code,
and the reasons live only in comments beside them: the `+1e-10` ridge
(`inst/julia/.../adjoint_ekf.jl:42-47`), the double ridge on `Pr` and `S`
(`kalman_filters.jl:466-486`), the Joseph form applied to the unridged
covariance, the dynamic-block restriction matching Stan's `derrind` subset
(`discrete_time_form.jl:156-158`), and `_laplace_popchol`'s reproduction of
`ctsm.stan:485-500` term for term. When the Stan path goes, those comments
become the only record of why the arithmetic is shaped this way, and the reverse
pass will be wrong if the forward's ridges are later tidied without changing
`_CTSEM_RIDGE`. Turn each into a plain statement of the convention before the
Stan source is removed, not after.

One deliberate divergence to record while both paths exist: the julia path caps
the unconstrained random-effect correlation coordinate at 5.2933 and Stan does
not, so on a design that cannot identify such a correlation the two backends
estimate over different parameter spaces.

## Sequence, when the time comes

1. Soft-deprecate: `backend='stan'` warns once per session, pointing at the
   julia setup. Names in the `ctStan*` family already alias the `ct*` names.
2. Stop compiling: drop the stan models from the build, keep the R-side stan
   fit *reading* code so old fit objects still summarise and plot.
3. Remove the stan fitting path and its dependencies; keep `ctStanFit` as a
   class name that old objects carry, or provide a converter.
4. Retire the reading code one release later.
