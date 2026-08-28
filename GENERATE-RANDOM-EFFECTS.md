# Generating individual differences and TI-predictor effects

Branch: `generateRandomEffects`, off `laplaceMultilevel`.

`ctGenerate(backend = "julia")` currently **refuses** a model with individually
varying parameters rather than silently generating a fixed-effects dataset from
it. That refusal is in `.ctGenerateJulia()` (`R/ctGenerateJulia.R`) and is
deliberate: a dataset with no between-subject variation, returned from a model
that asked for it, is a wrong answer dressed as a working one.

Rebase onto `laplaceMultilevel` before starting — the julia generation path
landed there and this worktree predates it.

## What the specification carries today

From `ctModel(n.TIpred = 2)$pars`:

| column | type | meaning |
|---|---|---|
| `indvarying` | logical | *whether* this parameter varies over subjects |
| `sdscale` | numeric (1) | multiplier on the population-SD **prior**, not an SD |
| `TI1_effect`, `TI2_effect` | logical | *whether* TI predictor *k* loads on this parameter |

and at model level: `tipredeffectscale`, `tipredsimputedscale`, `rawpopsdbase`,
`rawpopsdbaselowerbound`, `rawpopsdtransform`.

**Everything here is structure, and nothing is a value.** That is exactly right
for fitting — the population SD and each TI effect are themselves free
parameters with priors, so the fit only ever needs to know which ones exist.
Generation is the first consumer that needs the numbers, and there is no slot to
put them in. That, not the julia path, is the actual gap.

## The decision, made

**Values in the specification, priors as the fallback.** New per-row columns
carry a population SD and a TI effect size; anything left unset is drawn from
the prior that the fit would have used.

This was chosen over the two alternatives deliberately:

- *Priors only* (no spec change) would make `ctGenerate` a true prior-predictive
  sampler and cost nothing to build, but a user could not say "give me a
  population SD of 0.3" — which is the first thing anyone doing a power analysis
  wants, and the main reason to generate data at all.
- *Values only, error if unset* would be the most explicit, but `ctModel` sets
  `indvarying = TRUE` by default on T0MEANS and MANIFESTMEANS, so nearly every
  model would error until the user supplied numbers they may not care about.

The fallback is what keeps the default path working; the columns are what make
it controllable. Expect the columns to be the expensive half: `pars` is read by
`ctModel`, every matrix writer, both backends' preparation code, and anything
that prints or subsets it.

## Work, in order

1. **Add the columns.** A population SD per `indvarying` row, and a TI effect
   size per (row, predictor) pair. Follow the naming already there —
   `TI<k>_effect` is the logical, so the value wants a name that cannot be
   confused with it. Default `NA`, meaning "draw from the prior", which keeps
   every existing model behaving as it does now.
2. **Teach `ctModel()` to accept them** without breaking the existing calling
   conventions, and make sure a round trip through `ctModel` preserves them.
3. **Extend `.ctGenerateDefaults()`** in the same spirit as the fixed-effects
   values already there: chosen so generated data looks like data, documented as
   simulation defaults and not estimates.
4. **Draw the subject-level parameters** and apply the TI design, then hand the
   engine a per-subject parameter vector.
   `.ctJuliaPrepare(..., intoverpop = "augmented")` already builds the
   individual-effect structure the fit uses; generation needs the same structure
   populated with drawn values rather than estimated ones.
5. **Populate the TI predictor columns** in `.ctGenerateSkeleton()`, which
   currently writes zeros for them (and for TD predictors) — honest for a value
   the caller did not supply, useless once TI effects exist. `ctGenerate`'s own
   path uses `TDPREDMEANS`; there is no TI equivalent yet.
6. **Remove the refusal** in `.ctGenerateJulia()` and the `indvarying` workaround
   in the tests.
7. **Check the R path too.** `ctGenerate(backend = "r")` does generate individual
   differences; whether the two paths agree given a seed is worth knowing before
   claiming they are interchangeable.

## Watch for

- `ctModel` sets `indvarying = TRUE` by default on T0MEANS and MANIFESTMEANS
  rows, so the current refusal fires on models nobody thought were multilevel.
  Any partial implementation should narrow the refusal rather than widen it.
- The standard normals are drawn on the R side so `set.seed()` behaves. Subject
  parameter draws must come from the same stream, or seeding becomes a
  half-truth.
- Run tests with `Sys.setenv(NOT_CRAN = "true")`. The julia test files open with
  `skip_on_cran()`, and a bare `Rscript` run reports every skipped file as
  passing.
- `devtools::load_all(compile = FALSE)` for R-side iteration; a scratch install
  is only needed when C++ changes.
