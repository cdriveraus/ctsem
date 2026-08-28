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

## Status: implemented

Done, on this branch:

1. **`indvaryingsd`** on `pars`, and **`<TI>_effectsize`** beside each
   `<TI>_effect`. Both default to `NA` meaning unstated, so every existing model
   fits and generates exactly as before.
2. **The refusal is gone.** A model with varying parameters generates.
3. **Values are stated on the parameter's natural scale**, converted by the
   derivative of its own transform at the population mean. Exact for a mean
   parameter (`10 * param`, derivative constant); the delta method otherwise,
   and generation says which case it is in.
4. **TI predictors are drawn** — one standard normal per subject — so an effect
   has something to act on.
5. Unstated values use the prior's **centre**, not a draw from it, so two
   generations from one specification agree.

Measured: population sd 0.5 comes back as 0.493 and 2 as 1.94 on 400 subjects;
TI effects of 0.5 and 1.5 come back as slopes of 0.473 and 1.473.

### The finding that made it small

There was no sampler to write. Under `intoverpop = 'augmented'` a varying
parameter *becomes a state* — `mm` appears as row 2 of T0MEANS, its population
sd as `julia_popcov_2_2` in T0VAR — and `ctsem_generate` already draws the T0
state. Between-subject variation came out before any of this was written. The
work was a mapping, not a sampler.

### Three bugs it surfaced

- **Resolving free parameters destroyed the random effects.** Assigning a value
  makes a parameter fixed, and a fixed parameter is not augmented, so the
  population structure never existed: `npar = -Inf`, every stated spread
  ignored. Varying parameters stay free now.
- **`ctStanModelIntOverPop` built T0VAR rows positionally**, encoding the column
  order of `pars` as a literal. Any new column misaligned every field, and
  `rbind` warned rather than failed. Name-based now, and immune to the next
  column anyone adds.
- **`npar` was counted from the parameter table alone**, missing TI coefficients
  and the Laplace block, so the raw vector was short and the engine raised a
  `BoundsError` naming one of its own internals.

## Still to do

- The R and julia paths have not been compared for agreement given a seed.
  `ctGenerate(backend = 'r')` generates its own individual differences via
  `TRAITVAR` and `MANIFESTTRAITVAR`; whether the two produce the same spread
  from the same specification is worth knowing before calling them
  interchangeable.
- Correlations between random effects are in `spec$random_effects` with
  `type = 'correlation'` and are not settable — only the standard deviations
  are. A model with two varying parameters generates them uncorrelated.
- Nothing sets a population mean other than `.ctGenerateDefaults()`. A user
  wanting a varying parameter centred somewhere specific cannot say so, because
  the row must stay free and `value` is what "free" means.
- Multilevel (`indvarying_<idname>`) is untouched.

## Watch for

- `ctModel` sets `indvarying = TRUE` by default on T0MEANS and MANIFESTMEANS
  rows, so most models are multilevel without anyone asking.
- The standard normals are drawn on the R side so `set.seed()` behaves. Anything
  added must draw from the same stream.
- Run tests with `Sys.setenv(NOT_CRAN = "true")`. The julia test files open with
  `skip_on_cran()`, and a bare `Rscript` run reports every skipped file as
  passing.
- `devtools::load_all(compile = FALSE)` for R-side iteration. This worktree has
  no `src/`; copy `ctsem/src/*` into it once or `load_all` fails on the DLL.
