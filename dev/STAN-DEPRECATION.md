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

## Blockers

Things the julia path does not yet do that the stan path does. Each must be
resolved or explicitly dropped before deprecation starts.

| capability | stan location | status | found by |
|---|---|---|---|
| HMC via Stan (`optimize=FALSE`) | `R/ctFit.R` | julia has its own NUTS sampler; parity of diagnostics not yet reviewed | survey |
| `fullbootstrap` uncertainty | `R/ctOptimUncertainty.R` | not supported on julia per `JULIA-BACKEND.md` | survey |
| variational Bayes (`vb=TRUE`) | `R/ctFit.R` | refused on julia | survey |
| `summary(priorcheck=)` | `R/summary.ctStanFit.R` | accepted and ignored on julia | survey |

## Sequence, when the time comes

1. Soft-deprecate: `backend='stan'` warns once per session, pointing at the
   julia setup. Names in the `ctStan*` family already alias the `ct*` names.
2. Stop compiling: drop the stan models from the build, keep the R-side stan
   fit *reading* code so old fit objects still summarise and plot.
3. Remove the stan fitting path and its dependencies; keep `ctStanFit` as a
   class name that old objects carry, or provide a converter.
4. Retire the reading code one release later.
