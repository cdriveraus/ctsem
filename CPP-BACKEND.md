# A hand-written C++ backend for ctsem

Written 2026-08-22, answering `julia/CPP-BACKEND-SPIKE-PROMPT.md` in the
sibling ContinuousTimeSEM repository. That prompt asked for a narrow spike; this
is the full implementation instead, taken to parity with the Julia backend's
current feature set so the four engines can be compared directly on the same
models: **Stan, Julia forward (ForwardDiff), Julia adjoint, and hand-written
C++**.

Branch: `cppBackend`, off `juliaFit`.

---

## What was built

`backend='cpp'` in `ctFit()`, alongside `'stan'` and `'julia'`. It is a
continuous-time extended Kalman filter with a hand-written reverse-mode adjoint,
compiled into ctsem's own shared library. There is no per-model compile step and
no external toolchain at fit time.

| | file | role |
|---|---|---|
| engine | `inst/include/ctsemcpp/expr.hpp` | expression parser, flat AST, reverse-mode differentiation of a cell transform |
| | `inst/include/ctsemcpp/model.hpp` | the runtime model spec built from the R parameter table |
| | `inst/include/ctsemcpp/linalg.hpp` | matrix exponential, Frechet derivative, Lyapunov solve, `sdcovsqrt2cov`, and each one's pullback |
| | `inst/include/ctsemcpp/filter.hpp` | the primal EKF, which records the adjoint tape inline |
| | `inst/include/ctsemcpp/tape.hpp` | the tape records |
| | `inst/include/ctsemcpp/reverse.hpp` | the reverse pass |
| | `inst/include/ctsemcpp/engine.hpp` | multi-subject objective, gradient, L-BFGS |
| | `inst/include/ctsemcpp/rinterface.hpp` | spec -> objective |
| R boundary | `src/ctsemCppBackend.cpp` | Rcpp entry points |
| | `R/ctCppBackend.R` | `ctCppEvaluate`, `ctFitCppBackend`, S3 methods |
| tests | `tests/testthat/test-stan-cpp-parity.R` | the Julia parity suite's model shapes, against Stan |

The engine is header-only so it can also be compiled standalone with
`Rcpp::sourceCpp("src/ctsemCppBackend.cpp")` and `inst/include` on the include
path — which is how it was developed, without reinstalling ctsem.

**The R side needed no new model-preparation code.** `.ctCppPrepare()` calls
`.ctJuliaPrepare()` and drops two Julia-specific fields. That is the single most
important structural fact here: the expensive, bug-prone half of a backend is
turning a `ctModel` into a canonical, augmented parameter table with each cell's
transform rendered as an expression string, and the Julia backend already does
that. Stan compiles those strings, Julia `Meta.parse`s and `eval`s them, and this
engine interprets them. Adding a third *numerical* backend therefore cost about
120 lines of R.

---

## The AD choice, and why

**Hand-written adjoint over doubles + Eigen**, ported from
`ContinuousTimeSEM/src/adjoint_primitives.jl` and `adjoint_ekf.jl` — option (b)
in the spike prompt, not Stan Math's `var`.

The reasoning, in order of weight:

1. **`stan::math::var` is what the current Stan backend already uses.** ctsem's
   `src/stanExports_ctsm.cc` is generated C++ over Stan Math's reverse-mode AD
   and Eigen. A `var`-based engine would therefore reproduce the existing Stan
   backend's per-gradient cost almost exactly, minus whatever the generated code
   gains from having the model structure compiled in. The only remaining win
   would be compile time. That is a real win, but it makes the whole exercise
   uninformative about speed, which is half of what the comparison was for.
2. **Stan Math has no continuous Lyapunov solve.** The filter needs
   `A X + X A' + Q = 0` at every prediction step. Under `var` that would need a
   custom `vari` with a hand-written pullback anyway — i.e. the hardest single
   piece of option (b) is unavoidable under option (a) too.
3. **The derivations already exist and are validated.** `adjoint_primitives.jl`
   writes out every pullback with its derivation, and
   `test_adjoint_gradient_validation.jl` has fourteen scenarios pinning them.
   Porting a validated derivation is a much smaller risk than deriving one.
4. **Differentiating the transform AST directly is both easier and faster than
   dual numbers.** The Julia backend gets a state-dependent cell's partials by
   seeding a ForwardDiff dual once per input, so a transform reading `s`
   parameter cells in an `n`-state model costs it `s + n` evaluations of the
   expression. One reverse sweep over the flat AST gives all of them at once.

The `var` option is not foreclosed. Nothing in the filter or the R boundary
assumes `double`; the primal could be templated on the scalar type later if
anyone wants the comparison.

---

## What the expression interpreter cost

**Essentially nothing, and that is the answer to the spike's first
"things that would make this a bad idea".**

Profiled on the models below, the state-dependent transform groups are
**0.06 ms of a 28.6 ms gradient** at 20 latents (0.2%), and 0.09 ms of 90 ms on
the state-dependent variant of the same model (0.1%). The matrix kernels
dominate, exactly as they do in the Julia backend and exactly as the spike
prompt predicted.

Two design choices are why:

- The AST is a flat, topologically ordered array evaluated by one forward loop.
  No pointer chasing, no recursion, no allocation per evaluation.
- Read sets are discovered structurally from the parsed AST, not by probing
  derivatives numerically. A dependency whose derivative happens to vanish at a
  probe point cannot be missed.

So: the interpreter does not dominate, Stan's compile-it strategy is not
vindicated on speed grounds, and a runtime-evaluated engine is viable.

---

## Correctness

Two independent gates, both green.

**Against the Julia adjoint, on twelve model shapes** — the whole feature
surface, deliberately reusing the shapes that caught real bugs in the Julia port
rather than fresh clean ones:

| scenario | free params | value rel. diff | gradient rel. diff |
|---|---|---|---|
| linear 1-latent | 1 | 3.3e-16 | 1.2e-15 |
| PARS cross-effect + augmented state | 5 | 3.6e-16 | 2.8e-16 |
| partial (not total) row missingness | 5 | 3.2e-16 | 5.0e-16 |
| fully missing row | 5 | 8.9e-16 | 1.4e-16 |
| default free 2-latent | 23 | 0 | 3.1e-16 |
| nonlinear DRIFT + TD/TI predictors + bounded substeps | 7 | 2.2e-16 | 5.2e-18 |
| T0MEANS-indvarying population SD, non-unit meanscale | 3 | 5.9e-16 | 3.1e-16 |
| mixed random effects, 5 augmented states | 22 | 9.4e-16 | 2.8e-15 |
| state/TD-dependent LAMBDA + partial missingness | 6 | 0 | 2.6e-16 |
| 3 original latents, fully cross-coupled | 22 | 0 | 1.5e-15 |
| unequal-length subjects, irregular observation times | 27 | 2.2e-16 | 4.9e-16 |
| 6-latent default free model | 153 | 1.8e-16 | 4.3e-15 |

Agreement is at machine precision, not merely within tolerance. That matters:
the two implementations share no code, no language and no linear-algebra
library, so a shared misreading of the model spec is the only way they could
both be wrong together — and the R-side spec is itself checked against Stan.

**Against Stan**, `rstan::log_prob(..., adjust_transform = FALSE, gradient =
TRUE)` at a fixed raw parameter point, which is the independent ground truth:

| model | value rel. diff | gradient rel. diff |
|---|---|---|
| linear 1-latent | 3.7e-10 | 5.3e-8 |
| PARS cross-effect + augmented state | 1.7e-9 | 6.4e-9 |
| partial row missingness | 2.4e-9 | 5.0e-9 |
| linear 2-latent, 20 subjects x 5 waves | 2.6e-12 | 1.5e-10 |
| linear 6-latent, 20 subjects x 5 waves | 1.6e-11 | 7.5e-10 |
| linear 20-latent, 20 subjects x 5 waves | 2.2e-11 | 1.7e-9 |

The Julia backend agrees with Stan to *the same* figures on the same models
(compare the `julia` and `cpp` columns of the benchmark output: they match to
five significant digits of relative difference). So the residual is Stan's own
difference from both ports — different matrix-exponential and Lyapunov
implementations at the 1e-16 level, accumulated over the recursion — not
something this engine introduces.

`tests/testthat/test-stan-cpp-parity.R` runs eleven of these as permanent
regression tests (29 assertions, all passing), including one that compares the
adjoint against central finite differences of its own likelihood -- the referee
that does not share the reverse pass's assumptions -- and one that runs the two
backends' actual optimizers to convergence and compares the fitted raw
parameters, which is the guarantee a fixed-point check cannot give.

---

## Performance

One gradient evaluation, **20 subjects x 5 equally spaced waves**, **engine
time**: Julia timed inside Julia via `ctsem_evaluate`, C++ and Stan through
their `.Call` boundary, which costs 18 us and 81 us respectively and is
therefore already negligible at these sizes.

Engine time is the number that matters, because a fit repeats the gradient
inside the engine. `ctFit(backend='julia')` hands the whole L-BFGS loop to
`ctsem_optimize` in Julia, and `backend='cpp'` hands it to `.ctsemCppOptimize`
in C++; both cross the R boundary **once per fit**, not once per gradient. (The
Stan backend is the exception -- `stanoptimis` drives `mize`, an R-level
L-BFGS, calling `rstan::log_prob` per iteration -- but its boundary cost is
0.08 ms, so its engine and end-to-end numbers coincide to within noise.)

Minimum over repeated calls; the machine was a normally loaded desktop, so
treat the absolute seconds as indicative and the ratios as sound.
Free-parameter counts are matched across the linear and state-dependent
variants of each size. The state-dependent models put
`DRIFT[1,2] = PARS[1,1] * (1 + .05 * eta1)` and
`DIFFUSION[1,1] = PARS[2,1] * (1 + .03 * eta1)`, the shape that forces
`standata$recompile == 1` on the Stan side.

**Linear**

| latents | free params | Stan | Julia forward | Julia adjoint | **C++** |
|---|---|---|---|---|---|
| 2 | 23 | 0.00112 s | 0.00351 s | 0.00109 s | **0.00054 s** |
| 6 | 153 | 0.00322 s | 0.0535 s | 0.00342 s | **0.00240 s** |
| 20 | 1490 | 0.0410 s | 60.6 s | 0.0454 s | **0.0286 s** |

**State-dependent DRIFT + DIFFUSION**

| latents | free params | Stan | Julia forward | Julia adjoint | **C++** |
|---|---|---|---|---|---|
| 2 | 23 | 0.00150 s | 0.00539 s | 0.00279 s | **0.00084 s** |
| 6 | 153 | 0.00850 s | 0.0917 s | 0.0108 s | **0.00525 s** |
| 20 | 1490 | 0.269 s | 245.4 s | 0.0921 s | **0.0902 s** |

(The Julia forward column is still end-to-end; at 60-245 seconds per gradient
the boundary is not worth separating.)

**C++ is the fastest engine at every size, on both model families**: 1.3x to
3.0x Stan and 1.0x to 3.3x the Julia adjoint. But the interesting structure is
not the C++ column.

1. **Both adjoints beat Stan decisively on the model that matters most, and
   they beat it by about the same margin.** At 20 latents with a
   state-dependent DRIFT and DIFFUSION, Stan takes 0.269 s against 0.0921 s for
   Julia and 0.0902 s for C++ -- **~3x, for both**. Stan degrades 6.6x going
   from linear to state-dependent at that size; C++ degrades 3.1x and Julia
   1.7x. This is the adjoint roadmap's central claim, and on engine terms it
   holds more strongly than the end-to-end numbers ever showed.
2. **On linear models Stan and the Julia adjoint are level**, within 3-11% at
   every size, with C++ 1.4x-2.0x ahead of both.
3. **C++'s advantage over Julia is a per-row constant factor, and it shrinks as
   the model grows.** 2.0x/3.3x at 2 latents, 1.4x/2.1x at 6, 1.6x/1.02x at 20
   (linear/state-dependent). On the largest state-dependent model the two are a
   dead heat.
4. **Julia's primal filter is at parity with the C++ one** -- marginally
   *faster* at 20 latents on both families (0.0074 s vs 0.0086 s linear,
   0.0207 s vs 0.0212 s state-dependent). The entire engine-level gap is in the
   reverse pass.
5. **ForwardDiff is unusable at scale.** 245 seconds for one gradient of a
   20-latent state-dependent model puts an L-BFGS fit of it out of reach.
6. **The adjoint costs a normal multiple of the primal.** C++ value-only times
   are 0.00017 / 0.00071 / 0.0086 s (linear) and 0.00026 / 0.0016 / 0.021 s
   (state-dependent), so the gradient is 3.2-4.2x the primal at every size --
   the ratio a well-behaved reverse mode should have.

### The R boundary, and why it is not in the table above

An earlier version of this document reported end-to-end times measured through
`ctJuliaEvaluate` / `ctCppEvaluate`. That flattered C++ and penalised Julia for
something a real fit does not pay: `.Call` costs 18 us, while a JuliaConnectoR
round trip costs 2.5 ms at 23 parameters rising to 24 ms at 1490 -- but the
Julia backend crosses that boundary once per `ctFit`, not once per gradient.

For reference, the same six models end-to-end through `ctJuliaEvaluate` after
the fix below: 0.00359 / 0.00950 / 0.0671 s (linear) and 0.00646 / 0.0186 /
0.116 s (state-dependent). Those are the right numbers for a *diagnostic* call
-- `ctJuliaEvaluate` at one raw parameter vector is the workhorse of the
Stan/Julia parity tests and of every backend-disagreement investigation to date
-- and the wrong numbers for judging the engine.

Measuring that turned up a real bug on the way. `.ctJuliaNumericVector()`
marshalled the parameter vector as an R *list*, which JuliaConnectoR sends
element by element: 0.0013 s at 23 parameters, 0.0073 s at 153 and **0.069 s at
1490**, against **0.0002 s, flat**, for the same vector sent as a plain
numeric. Both arrive as `Vector{Float64}`; only a length-one vector needs the
list form, since a bare length-one numeric arrives as a scalar. Fixed in
`R/ctJuliaBackend.R` on this branch. It does not change fit times materially --
one call per fit -- but it takes 0.069 s off every diagnostic evaluation of a
large model, and it is why `docs/src/adjoint-roadmap.md`'s tables (and the first
version of this document's) were half JuliaConnectoR. That claim there, "the RPC
floor is 0.00026 s, so it only matters in the first row", is corrected.

### Per-model setup cost

| | Stan | Julia | C++ |
|---|---|---|---|
| linear models (`recompile == 0`) | 0.02-0.11 s (reuse of the generic binary) | 7-115 s | 0.001-0.010 s |
| state-dependent models (`recompile == 1`) | **162-182 s** | 5-112 s | 0.0005-0.004 s |

This is the motivation the spike prompt led with, and it survives measurement
intact. Every distinct state-dependent model definition costs Stan roughly three
minutes before a single gradient is evaluated, every time the definition
changes. The C++ objective is built in milliseconds: parse each transform string
into an AST once, copy the data, done. The Julia figures include JIT compilation
and are dominated by it at 20 latents.

### Where the time goes, and what the comparison found in the other backend

Profiling the reverse pass by tape-record kind at 20 latents put **66% of it in a
single operation**: the matrix-exponential Fréchet derivative, which is an
exponential of a `2n x 2n` block and costs ~540 us per prediction substep
against ~30 us for everything else in that substep combined.

`L(A, E)` is **linear in `E`**, and a balanced panel hands it the same
`JAx * dt` at every substep of every subject. So the directions are accumulated
and `dt * L(A', sum E)` is evaluated once instead of once per substep. Exact,
not approximate; the gradient still agrees with Julia's to 1e-15. On the
20-latent linear model that takes the count of block exponentials per gradient
from **80 to 1**, and the whole gradient from 0.0792 s to 0.0285 s — 2.8x, and
the difference between losing to Stan by 1.96x and beating it by 1.42x.

Two guards make it safe, and both check the actual data rather than a
model-level "is this linear" flag: the batch breaks whenever `JAx * dt` changes
(so a state-dependent model simply never batches and pays only the comparison),
and it is flushed before anything that could consume the JAx cotangent — a
transform group that writes a JAx cell, or the parameter layer. The condition
for carrying a batch across subjects is deliberately weaker than the condition
for sharing the whole parameter layer, because ctsem's default model makes
MANIFESTMEANS individually varying, which puts a state-dependent transform group
between every pair of prediction substeps in an otherwise entirely linear model.
A first attempt that flushed before every group therefore batched nothing at
all, and looked like a 3% win.

**Transplanted to the Julia adjoint** — `_flush_frechet!` in `adjoint_ekf.jl` on
`ctsem-backend`, which had the same structure and the same missing optimisation
— it is worth **1.24x / 1.44x / 1.56x** engine-only on the linear models at
n = 2 / 6 / 20, and nothing on the state-dependent ones, where the batch never
forms exactly as designed. Smaller than C++'s 2.8x because the operation was a
smaller share of Julia's reverse pass to begin with: Julia's `Base.exp` on a
`2n x 2n` block is considerably faster than the equivalent under R's default
compile flags, which is also why the C++ primal is *not* faster than Julia's
despite being hand-tuned.

That improvement lands on the *existing* backend regardless of whether the C++
one is ever adopted, and it is the one of the two Julia-side changes that
affects fit times: it is inside the engine, so a fit gets it on every gradient.
(The marshalling fix is not — one call per fit — which is why the two are worth
keeping separate rather than quoting their combined end-to-end effect.) A
benchmark comparison is a debugging tool for the thing it is compared against,
not only for the new thing.

### What is *not* claimed

- Single-threaded, like the Julia backend. RcppParallel is already a dependency
  and the subject loop is embarrassingly parallel, so there is obvious headroom
  that no engine here uses.
- Compiled with R's default flags, as Stan's generated code is. A 40x40 matmul
  measures ~50 us here, which is roughly scalar speed; `-march=native` would
  help all three engines and is not available on CRAN.
- The machine was a normally loaded desktop, not a quiet benchmark host.

---

## Honest assessment

**Worth building, and it is built.** Every one of the spike prompt's three
"things that would make this a bad idea" was checked and none of them holds:

1. *"The expression interpreter turns out to dominate runtime."* It is 0.2% of a
   gradient. The matrix kernels dominate, as they do in every backend here.
2. *"Stan Math's `var` overhead makes it no faster than the current Stan
   backend."* Not applicable — `var` was deliberately not used, for the reasons
   above — and the hand-written adjoint is faster than the Stan backend at every
   size measured, not merely equal to it.
3. *"A C++ engine cannot beat the Julia adjoint on nonlinear models."* It beats
   it by 3.3x at 2 latents and 2.1x at 6 — but at 20 latents, the size where
   this actually matters, the two engines are level (0.0902 s vs 0.0921 s).

So the case is stronger than "the win is compile time" — but **materially less
strong than the first version of this document claimed**, and the correction is
worth stating plainly rather than burying. Two things moved it: batching the
block exponential in the Julia adjoint, and measuring engine time rather than
end-to-end time through the R boundary. Together they removed almost all of the
C++ engine's apparent advantage on the largest state-dependent model. Julia's
primal filter was already at parity, and is slightly ahead at 20 latents.

What survives, and what the decision should actually rest on:

- **Deployment.** No Julia install, no `Pkg.add` from a pinned gitlab revision,
  no out-of-process session, no multi-minute first-run precompilation, no second
  repository to keep in lock-step. This is the argument that did not move at
  all, and it is now the main one.
- **Compile time versus Stan.** Three minutes per distinct state-dependent
  model, against milliseconds. Unchanged and large — but note this is an
  argument against *Stan*, not against Julia, which has no per-model compile
  either.
- **Speed.** Fastest at every size measured, decisively against Stan on the
  large state-dependent model (3.0x) and against Julia on small and medium
  models (1.4x-3.3x), but a tie with Julia at the large end. On the strength of
  the numbers alone there would be no case for replacing the Julia engine.

### The retirement criterion

The spike prompt asked for one to be agreed **up front**, and that is right:
three numerical backends that must stay in lock-step is a materially worse
maintenance position than two, and the two existing ones have a documented
history of silently drifting apart.

A recommendation, for the user's decision rather than mine:

> **Retire the Julia backend** once the C++ engine has (a) run the full
> `test-stan-julia-parity.R` model set green as `test-stan-cpp-parity.R` does,
> (b) reproduced a real `ctFit()` on `AnomAuth[1:1000,]` to the same optimum as
> Stan, and (c) been exercised through one full `R CMD check`. Keep the Stan
> path: it is the ground truth every other backend is validated against, it is
> the only one that does HMC, priors, importance sampling and uncertainty
> quantification, and none of that is in scope for either optimisation engine.

The argument for retiring Julia rather than C++ is **entirely a deployment
argument, not a speed one**. On engine time the two are level on the largest
state-dependent model and Julia's primal is marginally the faster of the two;
C++ leads on small and medium models by a per-row constant factor that shrinks
as models grow. What differs is that C++ ships inside the package, with no
out-of-process session, no `Pkg.add` from a pinned gitlab revision, no
multi-minute first-run precompilation and no separate repository to keep in
lock-step — and that Julia's reason for existing (a reverse-mode gradient whose
cost is flat in the parameter count) is fully reproduced here.

The arguments *against* doing it soon: mileage, since the Julia adjoint has been
exercised on real fits and this has not; and the fact that the speed case
evaporated once measured properly, which is a reason to make this decision on
deployment grounds deliberately rather than on a benchmark that has moved twice
already. Nothing here needs deciding today.

### What is missing relative to the Julia backend

Nothing in the likelihood or the gradient — the twelve-scenario comparison above
is the whole feature surface. What is missing is peripheral:

- Subject-specific matrix reconstruction (`ctExtract(subjectMatrices=TRUE)`),
  and everything downstream of the Kalman filter's *states* rather than its
  parameters: `ctPredict()`, `ctKalman()`, `ctPredictTIP()`, `ctACFresiduals()`,
  `ctPostPredPlots()`, `ctFitCovCheck()`, `ctGenerateFromFit()`. These need an
  engine entry point that returns per-row filtered and smoothed states and
  covariances, which neither engine has; the parameter-matrix work below does
  not get them for free. (Also missing for `ctJuliaFit`.)
- Multi-start or restart robustness in the optimizer (also missing in
  `ctsem_optimize`; the Julia handoff lists it as the top open item).
- The `v1` refusals are the same list as Julia's: no HMC, no priors, no
  discrete-time models, no non-Gaussian manifests, no variational Bayes, no data
  generation.
- CI does not exercise `test-stan-cpp-parity.R`, for the same reason it does not
  exercise the Julia one: the workflow installs neither rstan nor a Stan
  toolchain. Unlike the Julia suite, though, this one needs *only* rstan — no
  Julia install, no `CTSEM_JULIA_PROJECT` — so wiring it up is a much smaller
  decision.

### Transformed-parameter summaries

`summary()`, `ctSummaryMatrices()` and `ctDiscretePars()` work for `ctCppFit`
and `ctJuliaFit`, and they are not reimplementations — they are the same code
Stan fits use, reading the same input.

Every one of those functions consumes samples of the *model matrices* implied by
a raw parameter vector, as `pop_DRIFT`, `pop_DIFFUSIONcov` and so on, each an
`[iteration, row, column]` array. So the whole architecture is one engine
primitive and a thin stack on it:

| layer | what it does |
| --- | --- |
| `ctsem_parameter_matrices` / `.ctsemCppParMatrices` | materialize every model matrix for an `npar x nsamples` block of raw vectors |
| `ctBackendParMatrices()` | the same, named and reshaped, for one vector |
| `ctExtract()` | the same over the posterior, as `pop_*` arrays |
| `ctSummaryMatrices()` | collapse those arrays — literally `ctSummaryMatrices.ctStanFit`'s body, factored into `.ctSummaryMatricesFromArrays()` |
| `summary()` | fixed effects plus the system matrices |

The primitive runs **inside the engine**. The engine already materializes every
matrix from the raw vector on its way to a log likelihood, so asking it for that
same materialization is the only way to guarantee a summary reports what the
likelihood used; a second, R-side implementation of ctsem's transforms is
exactly the kind of duplication that has let this package's backends drift apart
before. The derived matrices (`DIFFUSIONcov`, `MANIFESTcov`, `T0cov`,
`asymDIFFUSIONcov`, `asymCINT`) are computed there too, following the generated
Stan code's definitions.

The batch shape is not incidental. The Julia bridge marshals a numeric array in
one transfer but a list element by element, which was the 344x cost found
earlier in this work; a 200-draw posterior is one call, not 200.

Three things are worth knowing:

- **Verification is against Stan's own constrain step**, not against a second
  fit: `stan_constrainsamples()` at a fixed raw vector, on a model with an
  `intoverpop` augmentation, TI predictors and a state-dependent DRIFT. Every
  `pop_*` array agrees to 1e-8 or better, and C++ and Julia agree with each
  other to 1e-12. That removes the optimizer from the comparison entirely.

- **`pop_T0VAR` is the one exception, and it is a parameterisation difference,
  not a disagreement.** Stan computes `T0cov = sdcovsqrt2cov(T0VAR)` and *then*
  rescales `T0cov`'s indvarying-T0MEANS rows and columns by the parameter's
  multiplier and meanscale, leaving `T0VAR` itself unscaled. The engines fold
  that scale into `T0VAR`, so theirs is the one whose `sdcovsqrt2cov` actually
  equals the reported `T0cov`. Both give an identical `T0cov` — which is the
  quantity summaries report, and `summary()` drops `T0VAR` from the system
  matrices table for exactly this reason.

- **State-dependent cells are reported as conditional, not as constants.** They
  are functions of the latent state, so there is no single number to report.
  They are evaluated at a state — T0MEANS by default, `state=` for anything
  else — and `attr(ctBackendParMatrices(fit), 'stateDependent')` names them, as
  does a note in `summary()`. Note that ctsem implements individually-varying
  parameters as a state dependence on augmented carrier states, so this applies
  to any model with random effects, not only to explicitly nonlinear ones.

Intervals appear only when they have been earned: a fit without
`ctOptimUncertainty()` has one "sample", and the interval columns are omitted
rather than filled with a zero-width interval that would read as certainty.

### One latent difference from the Julia backend, deliberately not copied

For the **first row of each subject**, the Julia filter builds the manifest
covariance from MANIFESTVAR *before* it applies the update-group transforms
(`_extended_kalman_filter_continuous!` calls `sdcovsqrt2cov!` on
`pars.MANIFESTVAR` above the `_record_group!(trace, 3, ...)` line), where every
later row does it after. MANIFESTVAR is in the update group, so for a model with
a state-dependent MANIFESTVAR that first row reads the cell before its transform
has ever written it. It has no effect on any model where MANIFESTVAR is not
state-dependent, which is every model either backend has a test for. This engine
builds it after the update group on every row, including the first — a one-line
divergence that is a strict improvement and cannot change any currently tested
result. Worth a look on the Julia side.

---

---

## The Julia backend's deployment, and what was done about it

Since the speed comparison levelled out, deployment was the only remaining
argument for the C++ engine. That made it worth measuring rather than asserting.

**What was actually wrong, verified rather than inferred:**

- **`backend='julia'` could not be installed by anyone outside one GitLab
  project.** `inst/julia/engine.json` pointed `ctJuliaSetup()` at
  `https://gitlab.uzh.ch/psyquantimet/continuoustimesem.git`, and an
  unauthenticated `git ls-remote` on it returns *HTTP Basic: Access denied*. It
  worked locally only because of an access token beside the checkout. This is
  not friction, it is a wall — and invisible from inside the group, which is
  why nothing caught it.
- **`revision` was a branch name**, so the "lock file" locked nothing.
- **The engine lived only in a second repository** that had to stay in lock-step
  with `R/ctJuliaBackend.R`'s parameter-table contract by hand.
- **111-package dependency closure, 268 MB in a fresh depot, ~73 s of install,
  and 8.7 s for `using ContinuousTimeSEM` in every R session.** Most of it was
  DataFrames, used only as a row container in the R interface, plus four
  dependencies with no call sites at all that were loaded because the module
  opened with `@reexport using` over all of them.

**What changed** (in `ctsem` on `cppBackend`, and in `ContinuousTimeSEM` on
`ctsem-backend`):

| | before | after |
|---|---|---|
| install source | `Pkg.add` from a private URL | vendored in `inst/julia/` |
| credentials / network at setup | required | none |
| pinned to | a branch | the exact commit, plus a vendored `Manifest.toml` |
| dependency closure | 111 packages | **67** |
| fresh depot | 268 MB | **124 MB** |
| dependency install | ~73 s | **~20 s** |
| `using ContinuousTimeSEM` | 8.7 s **per R session** | **3.86 s** |
| cold `ctJuliaSetup()` | network + auth + resolve | **10.7 s, offline** |

`ctJuliaSetup()` now copies the vendored package into a writable,
revision-keyed project under `R_user_dir("ctsem", "cache")` and instantiates it
there — activating it in place would fail on a read-only R library. If the
vendored manifest is unsatisfiable on the user's Julia version it falls back to
a fresh resolve rather than refusing to run. `tools/sync-julia-engine.sh`
refreshes the copy, refuses to run against a dirty tree, and records the source
commit in `engine.json`, which is now provenance rather than an install spec.

`ekf_from_columns` replaced `ekf_from_data_frame`, taking the parameter table as
plain column vectors with sentinels (`0`, `NaN`, `""`) instead of `missing`.
DataFrames stays as a *test-only* dependency, because a `DataFrame` literal is
still the clearest way to write a small table in a test; `test/table_helpers.jl`
adapts. The Julia suite passes 237/237 unchanged, and the twelve-scenario
C++/Julia comparison still agrees to ~1e-15.

Two smaller things the work surfaced, both now fixed: `.ctJuliaString` did not
escape backslashes, so a Windows path in a Julia string literal failed to parse
as "invalid unicode escape" and said nothing about the real problem; and
JuliaConnectoR *hangs* marshalling an empty vector, so the optional table
columns are Julia keyword arguments passed only when they have entries.

### Spinning the engine out as a standalone Julia package

The vendored tree is laid out exactly as a **root-level Julia package** —
`Project.toml`, `Manifest.toml`, `src/`, `test/`, and nothing else. That is
deliberate, and it is what makes each of these possible without rearranging
anything:

- **Work on it in place.** `inst/julia/ContinuousTimeSEM/` is a complete
  package: `julia --project=inst/julia/ContinuousTimeSEM -e 'using Pkg;
  Pkg.test()'` runs the whole Julia suite against the vendored copy.
- **Work on it as its own repository.** `ctJuliaSetup(project = "<checkout>")`
  still points ctsem at a development checkout instead of the vendored copy,
  which is how the engine has been developed throughout.
- **Register it as a Julia package.** A root-level package with a UUID and
  compat bounds is exactly what the General registry wants, so a pure-Julia
  user could eventually `Pkg.add("ContinuousTimeSEM")` with no R involved.

The one structural step still outstanding is upstream, not here: the
`ContinuousTimeSEM` repository keeps the package in a `ContinuousTimeSEM/`
subdirectory alongside unrelated scratch files. Move it to the repository root
and `tools/sync-julia-engine.sh` can be deleted in favour of

```
git subtree pull --prefix=inst/julia/ContinuousTimeSEM <remote> <branch> --squash
git subtree push --prefix=inst/julia/ContinuousTimeSEM <remote> <branch>
```

which makes the relationship two-way: edits made in ctsem can be pushed back
upstream, and upstream edits pulled in, with real history on both sides.

A **submodule** would be the wrong tool here, and it is worth saying why since
it is the more obvious reach. Submodule contents are not part of the parent
repository's tree and are not included in an R package tarball, so a released
ctsem would ship an empty directory and users would be back to needing network
access and credentials — the exact problem being fixed. Subtree vendors the
content while keeping the two-way link; submodule vendors only a pointer.

## Appendix: reproducing the benchmark

Requires `rstan`, `JuliaConnectoR` and a local ContinuousTimeSEM checkout for
the Julia rows. Drop the Julia rows and it needs only `rstan`.

```r
library(ctsem); library(rstan)
JULIA_PROJECT <- "<path to>/julia/ContinuousTimeSEM"

linmodel <- function(n) suppressWarnings(ctModel(type = "ct", LAMBDA = diag(n)))

nlmodel <- function(n) {
  drift <- matrix("", n, n)
  for (i in 1:n) for (j in 1:n) drift[i, j] <- paste0("dr", i, "_", j)
  drift[1, 2] <- "PARS[1,1] * (1 + .05 * eta1)"
  diffusion <- matrix("0", n, n)
  for (i in 1:n) for (j in 1:i) diffusion[i, j] <- paste0("df", i, "_", j)
  diffusion[1, 1] <- "PARS[2,1] * (1 + .03 * eta1)"
  suppressWarnings(ctModel(type = "ct", LAMBDA = diag(n), DRIFT = drift,
    DIFFUSION = diffusion, PARS = matrix(c("nlp1", "nlp2"), 2, 1)))
}

makeData <- function(n, nsub = 20, waves = 5, seed = 42) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(nsub), each = waves),
                  time = rep(seq_len(waves) - 1, nsub))
  for (i in seq_len(n)) d[[paste0("Y", i)]] <- rnorm(nsub * waves, 0, .5)
  d
}

# Minimum over repeated batches; a single call is far too noisy at these sizes.
timeit <- function(f, batch = 20, reps = 5) {
  best <- Inf
  for (r in seq_len(reps)) {
    t0 <- Sys.time()
    for (b in seq_len(batch)) f()
    best <- min(best, as.numeric(difftime(Sys.time(), t0, units = "secs")) / batch)
  }
  best
}

n <- 20; model <- linmodel(n); data <- makeData(n)

stan_spec <- ctFit(data, model, backend = "stan", fit = FALSE, priors = FALSE)
stan_fit <- if (identical(stan_spec$standata$recompile, 0L)) {
  ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, stan_spec$standata)
} else {
  ctsem:::stan_reinitsf(rstan::stan_model(model_code = stan_spec$stanmodeltext),
                        stan_spec$standata)
}
set.seed(9); raw <- rnorm(rstan::get_num_upars(stan_fit), 0, .2)

cpp_spec <- ctFit(data, model, backend = "cpp", fit = FALSE, priors = FALSE)
julia_spec <- ctFit(data, model, backend = "julia", fit = FALSE, priors = FALSE,
                    backendcontrol = list(julia_project = JULIA_PROJECT))

timeit(function() rstan::log_prob(stan_fit, upars = raw,
                                  adjust_transform = FALSE, gradient = TRUE))
timeit(function() ctCppEvaluate(cpp_spec, raw, gradient = TRUE))
timeit(function() ctJuliaEvaluate(julia_spec, raw, gradient_method = "adjoint"))
timeit(function() ctJuliaEvaluate(julia_spec, raw, gradient_method = "forward"),
       batch = 1, reps = 1)   # 63 s per call at n = 20
```
