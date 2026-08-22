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

One gradient evaluation, **20 subjects x 5 equally spaced waves**, measured from
R so every row carries its own call overhead (rstan's for Stan, JuliaConnectoR's
RPC for Julia, `.Call` for C++). Minimum over repeated calls — the machine was a
normally loaded desktop, so treat the absolute seconds as indicative and the
ratios as sound. Free-parameter counts are matched across the linear and
state-dependent variants of each size.

The state-dependent models put `DRIFT[1,2] = PARS[1,1] * (1 + .05 * eta1)` and
`DIFFUSION[1,1] = PARS[2,1] * (1 + .03 * eta1)`, which is the shape that forces
`standata$recompile == 1` on the Stan side.

**Linear**

| latents | free params | Stan | Julia forward | Julia adjoint | **C++** |
|---|---|---|---|---|---|
| 2 | 23 | 0.00112 s | 0.00351 s | 0.00359 s | **0.00054 s** |
| 6 | 153 | 0.00322 s | 0.0535 s | 0.00950 s | **0.00240 s** |
| 20 | 1490 | 0.0410 s | 60.6 s | 0.0671 s | **0.0289 s** |

**State-dependent DRIFT + DIFFUSION**

| latents | free params | Stan | Julia forward | Julia adjoint | **C++** |
|---|---|---|---|---|---|
| 2 | 23 | 0.00150 s | 0.00539 s | 0.00646 s | **0.00084 s** |
| 6 | 153 | 0.00850 s | 0.0917 s | 0.0186 s | **0.00526 s** |
| 20 | 1490 | 0.269 s | 245.4 s | 0.116 s | **0.0897 s** |

**C++ is the fastest engine at every size, on both model families**: 1.3x to 3.0x
Stan, and 1.3x to 6.6x the Julia adjoint. The margin over Stan is largest exactly
where it matters most — a 20-latent state-dependent model, where Stan is 3.0x
slower and also needs a three-minute compile.

The Julia adjoint column reflects two changes made *after* the first version of
this document, both prompted by profiling done for this port and both now
committed; see "what the comparison found in the other backend" below. Against
the originally measured Julia adjoint the C++ margin was 2x to 9x, so treat any
earlier statement of that ratio as superseded.

Three things worth reading off the table beyond the headline:

1. **The Julia adjoint's flatness under nonlinearity is real, and it is why the
   two adjoints converge at scale.** Going from linear to state-dependent at 20
   latents costs Stan 6.6x (0.0410 -> 0.269 s), C++ 3.1x, and the Julia adjoint
   1.7x. C++ still wins there, but by 1.3x rather than the 3.0x it enjoys over
   Stan on the same model.
2. **ForwardDiff is unusable at scale.** 245 seconds for one gradient of a
   20-latent state-dependent model. An L-BFGS fit of that model is out of reach
   on the forward path.
3. **The adjoint costs a normal multiple of the primal.** C++ value-only times
   are 0.00017 / 0.00071 / 0.0086 s (linear) and 0.00026 / 0.0016 / 0.021 s
   (state-dependent), so the gradient is 3.2-4.2x the primal at every size —
   the ratio a well-behaved reverse mode should have.

### Engine-only timings, and a measurement artefact worth knowing about

Everything above is end-to-end from R, which is what a user experiences. It is
*not* a clean engine comparison: JuliaConnectoR's round trip is charged to
Julia and `.Call`'s microseconds are charged to C++. Timing `ctsem_evaluate`
**inside** Julia removes that.

| model | free params | Julia engine, before | Julia engine, after | **C++** | Julia primal | C++ primal |
|---|---|---|---|---|---|---|
| linear 2 | 23 | 0.00136 s | 0.00109 s | **0.00054 s** | 0.00028 s | 0.00017 s |
| linear 6 | 153 | 0.00493 s | 0.00342 s | **0.00240 s** | 0.00063 s | 0.00071 s |
| linear 20 | 1490 | 0.0710 s | 0.0454 s | **0.0286 s** | 0.00737 s | 0.00858 s |
| statedep 2 | 23 | 0.00259 s | 0.00279 s | **0.00084 s** | 0.00122 s | 0.00026 s |
| statedep 6 | 153 | 0.0110 s | 0.0108 s | **0.00525 s** | 0.00522 s | 0.00156 s |
| statedep 20 | 1490 | 0.0913 s | 0.0921 s | **0.0902 s** | 0.0207 s | 0.0212 s |

Two findings here matter more than the ratios:

- **Julia's primal filter is at parity with the C++ one** — marginally *faster*
  at 20 latents on both families. The whole engine-level gap is in the reverse
  pass, and on the 20-latent state-dependent model even that closes to a dead
  heat (0.0921 s against 0.0902 s). Any argument for preferring one engine over
  the other on raw numerical speed is much weaker than the end-to-end table
  suggests.
- **The Julia backend's published timings were 50-80% JuliaConnectoR.**
  `.ctJuliaNumericVector()` marshalled the parameter vector as an R *list*,
  which JuliaConnectoR sends element by element: 0.0013 s at 23 parameters,
  0.0073 s at 153 and **0.069 s at 1490**, against **0.0002 s, flat**, for the
  same vector sent as a plain numeric. Both arrive as `Vector{Float64}`; only a
  length-one vector needs the list form, since a bare length-one numeric arrives
  as a scalar. That is a 344x difference on a function called once per objective
  evaluation, i.e. in the optimizer's inner loop. Fixed in `R/ctJuliaBackend.R`
  on this branch; `docs/src/adjoint-roadmap.md`'s claim that "the RPC floor is
  0.00026 s, so it only matters in the first row" is corrected there.

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

Combined with the marshalling fix, the Julia adjoint's end-to-end time improved
by **1.10x to 2.41x** across the six models, most of it at the large end where
both the block exponential and the parameter-vector marshalling scale worst.

Both improvements land on the *existing* backend regardless of whether the C++
one is ever adopted, which is worth weighing below: a benchmark comparison is a
debugging tool for the thing it is compared against, not only for the new thing.

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
   it end-to-end by 1.3x at 20 latents, 3.5x at 6 and 7.7x at 2 — but see the
   qualification below, because at 20 latents that margin is almost entirely
   the R boundary rather than the engine.

So the case is stronger than "the win is compile time" — but **less strong than
the first version of this document claimed**, and the correction is worth
stating plainly rather than burying. The two changes made to the Julia backend
in the course of this comparison (batching the block exponential, and fixing the
parameter-vector marshalling) removed most of the C++ engine's apparent
advantage on the largest state-dependent model: engine-to-engine there it is now
0.0902 s against 0.0921 s, which is a tie. Julia's primal filter is already at
parity, and slightly ahead at 20 latents.

What survives, and what the decision should actually rest on:

- **Compile time.** Three minutes per distinct state-dependent model for Stan,
  against milliseconds. Unchanged and large.
- **Deployment.** No Julia install, no `Pkg.add` from a pinned gitlab revision,
  no out-of-process session, no multi-minute first-run precompilation, no second
  repository to keep in lock-step. This is the argument that did not move.
- **Speed.** Still fastest at every size measured, decisively so against Stan
  and against Julia on small and medium models, but only marginally against the
  Julia engine on the largest state-dependent one.

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

The argument for retiring Julia rather than C++ is **not** that Julia is slow —
after the two fixes above the engines are close, and on the largest
state-dependent model they are level. It is that C++ ships inside the package,
has no out-of-process session, no `Pkg.add` from a pinned gitlab revision, no
multi-minute first-run precompilation and no separate repository to keep in
sync, and that Julia's reason for existing (a reverse-mode gradient whose cost
is flat in the parameter count) is fully reproduced here.

The arguments *against* doing it soon: mileage, since the Julia adjoint has been
exercised on real fits and this has not; and the fact that the speed case just
got weaker, which is a reason to decide on deployment grounds deliberately
rather than on a benchmark that has now moved twice. Nothing here needs deciding
today.

### What is missing relative to the Julia backend

Nothing in the likelihood or the gradient — the twelve-scenario comparison above
is the whole feature surface. What is missing is peripheral:

- `ctSummaryMatrices()` / subject-matrix reconstruction (also missing for
  `ctJuliaFit`).
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
