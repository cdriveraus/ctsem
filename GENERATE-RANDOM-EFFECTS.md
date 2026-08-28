# Generating individual differences and TI-predictor effects

Branch: `generateRandomEffects`, off `laplaceMultilevel`.

This branch exists because `ctGenerate(backend = "julia")` currently **refuses**
a model with individually varying parameters rather than silently generating a
fixed-effects dataset from it. That refusal is in
`.ctGenerateJulia()` (`R/ctGenerateJulia.R`) and is deliberate: a dataset with
no between-subject variation, returned from a model that asked for it, is a
wrong answer dressed as a working one.

Rebase onto `laplaceMultilevel` before starting -- the julia generation path
landed there in `437bfbc0` and this worktree predates it.

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
for fitting -- the population SD and each TI effect are themselves free
parameters with priors, so the fit only ever needs to know which ones exist.
Generation is the first consumer that needs the numbers, and there is no slot to
put them in. That, not the julia path, is the actual gap.

## The decision to make first

Two coherent answers, and they are not exclusive:

**(a) Draw from the priors.** Generation samples each population SD from
`rawpopsdbase`/`sdscale` and each TI effect from `tipredeffectscale`, using
exactly the priors the fit would use. No specification change; `ctGenerate`
becomes a genuine prior-predictive sampler, and generated-then-fitted recovery
becomes a real check of the prior. But a user cannot say "give me a population
SD of 0.3", which is the first thing anyone doing a power analysis wants.

**(b) Carry values in the specification.** New per-row columns (an
`indvaryingValue`, a TI effect value per predictor) alongside the existing
logicals. Explicit and controllable, at the cost of widening `pars` -- which
`ctModel`, every matrix writer, both backends' preparation code, and anything
that subsets or prints `pars` all touch.

A plausible resolution is (b) for control with (a) as the default fill, so an
unspecified model still generates something sensible and a specified one
generates what was asked. **Confirm with the user before building either** --
the choice determines how far the change reaches, and (b) is the one that
"will probably break lots of things".

## Then

1. Extend `.ctGenerateDefaults()` in the same spirit: values chosen so generated
   data looks like data, documented as simulation defaults and not estimates.
2. Draw subject-level parameters, apply the TI design, and hand the engine a
   per-subject parameter vector. `.ctJuliaPrepare(..., intoverpop = "augmented")`
   already builds the individual-effect structure the fit uses; the generation
   entry point needs the same structure populated with drawn values rather than
   estimated ones.
3. Populate the TI predictor columns in `.ctGenerateSkeleton()`, which currently
   writes zeros for them (and for TD predictors) -- honest for a value the
   caller did not supply, useless once TI effects exist. `ctGenerate`'s own path
   uses `TDPREDMEANS`; there is no TI equivalent yet.
4. Remove the refusal in `.ctGenerateJulia()` and the `indvarying` workaround in
   the tests.
5. Check the R path too. `ctGenerate(backend = "r")` does generate individual
   differences; whether the two paths agree given a seed is worth knowing before
   claiming they are interchangeable.

## Watch for

- `ctModel` sets `indvarying = TRUE` by default on T0MEANS and MANIFESTMEANS
  rows, so the refusal fires on models nobody thought were multilevel. Any
  partial implementation should narrow the refusal rather than widen it.
- The standard normals are drawn on the R side so `set.seed()` behaves. Subject
  parameter draws must come from the same stream, or seeding becomes a
  half-truth.
