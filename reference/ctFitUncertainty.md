# Compute or sample a fit's uncertainty

Recomputes the uncertainty of an optimised
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) fit and
redraws its approximate parameter samples, or samples the posterior by
Markov chain Monte Carlo. `ctFit` calls this function at the end of an
optimised fit. `ctOptimUncertainty` is the previous name, kept with the
arguments it had in ctsem 3.11.1.

## Usage

``` r
ctFitUncertainty(
  fit,
  method = c("hessian", "checked", "is", "sample", "sandwich", "opg", "bootstrap",
    "fullbootstrap", "stored"),
  intoverpop = NULL,
  filter = c("kalman", "particle"),
  checks = TRUE,
  finishsamples = NULL,
  cores = getOption("mc.cores", 2L),
  control = list(),
  verbose = 0,
  ...
)

ctOptimUncertainty(
  fit,
  uncertainty = c("hessian", "checked", "is", "bootstrap", "fullbootstrap", "sandwich",
    "opg", "stored"),
  draws = c("auto", "normal", "empirical", "imis"),
  finishsamples = NULL,
  cores = getOption("mc.cores", 2L),
  control = list(),
  verbose = 0,
  ...
)
```

## Arguments

- fit:

  An optimised `ctStanFit` or `ctJuliaFit`. `method = 'sample'` also
  accepts a julia fit that has been sampled.

- method:

  How the posterior is described: `'hessian'` (the default),
  `'checked'`, `'is'`, `'sample'`, `'sandwich'`, `'opg'`, `'bootstrap'`,
  `'fullbootstrap'` or `'stored'`. See the Methods and Sampling
  sections.

- intoverpop:

  How the random effects are integrated, with the values of
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)'s
  argument. `NULL` (the default) keeps the fit's own. `'quadrature'`,
  for a Laplace or variational fit, refines the fit by integrating
  numerically the subjects or groups the random-effects check flags, and
  describes the posterior there, as `ctFit(intoverpop = 'quadrature')`
  does. It works with `'hessian'`, `'checked'` and `'is'`. With
  `method = 'sample'`, `'none'` samples the random effects together with
  the parameters, and `'laplace'` samples the Laplace marginal (see
  Sampling).

- filter:

  How the latent states are integrated: `'kalman'` (the default), the
  fit's Kalman filter, or `'particle'`, a particle filter (see
  [`ctParticleLik`](https://cdriveraus.github.io/ctsem/reference/ctParticleLik.md)).
  `'particle'` works with `method = 'is'`, on fits whose random effects
  are not integrated by Laplace (see Methods).

- checks:

  Which checks to run, where the method does not replace them itself:
  `TRUE` (the default) for all, `FALSE` for none, or any of
  `'randomeffects'`, `'intervals'` and `'maximum'`. The random-effects
  check (julia, `intoverpop = 'laplace'`) scores the Laplace
  approximation for each subject or group, and gives a message
  suggesting `intoverpop = 'quadrature'` where integrating the random
  effects more exactly would move an estimate, or a warning where it
  would move one by a standard error or more. A fit that has been
  checked already is not checked again. With `method = 'sample'` on the
  joint posterior, a fit whose check finds a move of a standard error or
  more is refined first (see Sampling). The interval check (julia,
  `method = 'hessian'`) compares the Hessian's normal approximation with
  the posterior's gradients, records the result in
  `$uncertainty$details$normalcheck`, and reports in a message and in
  [`summary()`](https://rdrr.io/r/base/summary.html) when the posterior
  differs from the normal, suggesting `method = 'sample'`. Both checks
  usually take seconds. With `intoverpop = 'quadrature'` the interval
  check runs on the refined fit, and on nested models can take longer
  than sampling. `checks = 'randomeffects'` leaves it out. The maximum
  check (julia, `method = 'sample'` on the joint posterior) looks at
  where the sampler's first draws are centred. If that is far from the
  estimate, it optimises from there, and gives a warning if this reaches
  a higher log posterior than the estimate's. The point is kept in
  `$sample$maximumcheck$raw`. It costs nothing unless the draws are far
  from the estimate.

- finishsamples:

  Number of approximate raw parameter draws. If `NULL`, the number
  already in the fit's raw posterior (`fit$stanfit$rawposterior` for
  stan, `fit$estimate$rawposterior` for julia) is used, or 1000 if there
  is none. Not used by `'sample'`, which uses `control$draws`.

- cores:

  The most CPU cores the call uses at once, counting every process it
  starts. On a `ctJuliaFit` this sets the number of engine threads, up
  to the number the Julia session was started with, and with
  `method = 'sample'` and `processes = TRUE` the cores are shared
  between the worker processes. On a `ctStanFit`, more than one core
  starts R worker processes, which are used only for `'fullbootstrap'`
  and for `'is'` on data with a few thousand rows or more. Results are
  reproducible at `cores = 1`, and agree to rounding error above it.

- control:

  List of method specific settings.

  `accuracy` (julia, every method except `'sample'`) is a target
  accuracy in posterior sds. After the method runs, its own error
  estimate is compared with it (for `'checked'` and `'is'`, the Stein
  check's). If the method falls short, or has no error estimate
  (`'hessian'` and the score methods), the fit is sampled to that
  accuracy instead, with `minESS` and `rhatTarget` set from it: 0.5
  gives an effective sample size of about 125 and R-hat 1.07, and 0.3
  about 350 and 1.025. The method finally used is the one reported, and
  `$uncertainty$escalation` records what was done.

  For `method = 'sample'` the entries are the sampler settings in the
  Sampling section: `chains`, `warmup`, `draws`, `seed`, `saveEffects`,
  `processes`, `tipredMissing`, `sampler`, `particles`, `coreset`,
  `placement`, `minESS`, `meanESS`, `rhatTarget`, `maxDraws`,
  `settleTol`, `init_scale`, and for NUTS `maxdepth`, `target_accept`,
  `adapt_metric` and `adapt_effects`. The entries below do not apply to
  it.

  For the other methods, with their defaults: `ridge` (1e-8),
  `hessianStep` (1e-3), `bootstrapFitCores` (1), `bootstrapTol` (1e-5),
  `imisMaxIter` (50, or 100 for `'checked'`), `imisScaleInit` (1.1, or
  1.5 for `'is'` on julia), `imisTailScale` (1.1, or 1.2 for `'is'` on
  julia), `imisDf` (5, the degrees of freedom of the proposal's t
  components, `Inf` for normal components), `isESS` (200) and
  `isitersize` (1000, or 100 for `'checked'`). `imisBudget` limits the
  number of log density evaluations importance sampling uses (no limit
  by default). `imisComplement = TRUE` (`'checked'`) integrates over the
  other directions when weighting the sampled ones. This is more
  accurate where the spread of the other parameters changes with the
  sampled ones, as with the scale of a large covariance, and takes up to
  seven times as many evaluations. On julia, importance sampling
  evaluates a weighted subset of the subjects and corrects its weights
  where the subset is inaccurate, and `coreset = FALSE` evaluates every
  subject instead.

  When the check of `'checked'` finds departures from the normal, the
  normal is first moved to the posterior's mean and covariance by doubly
  adaptive importance sampling (van den Boom, Cremaschi and Thiery 2026;
  julia only) and checked again. `dais = FALSE` skips this step,
  `daisDraws` sets its draws per iteration (by default the number of
  parameters, between 64 and 256) and `daisMaxIter` its number of
  iterations (150).

  `parsteps` keeps stepwise fixed raw parameters fixed (stan only).
  `analyticHessian = FALSE` (julia only) uses a finite difference
  Hessian. `flat = 'natural'` (`'checked'` and `'is'`, julia only) gives
  parameters without a prior a flat prior on their natural scale, for
  example on a standard deviation itself where the default is flat on
  its log, so that a likelihood that levels off as an sd goes to zero
  still gives a proper posterior. The estimate does not change. `nodes`,
  `particles`, `substeps`, `transition`, `seed` and `correctEstimate`,
  for `intoverpop = 'quadrature'` and `filter = 'particle'`, are
  described under Methods. If the Hessian covariance needs numerical
  repair, a warning is given, with details in
  `$uncertainty$details$covariance`.

- verbose:

  Integer controlling progress detail.

- ...:

  Not used. Anything passed here is an error, since method settings go
  in `control`.

- uncertainty, draws:

  `ctOptimUncertainty()` only, with the argument names of ctsem 3.11.1.
  `uncertainty` is `method`, and `draws` (`'auto'`, `'normal'`,
  `'empirical'` or `'imis'`) sets the kind of draws the method makes.

## Value

The fit, of the class it came in as. For every method except `'sample'`,
the method, kind of draws, number of draws, cores and other settings
used are recorded in `fit$stanfit$uncertainty$settings` for a
`ctStanFit` and in `fit$uncertainty$settings` for a `ctJuliaFit`. See
Backend differences for the other parts each writes. For
`method = 'sample'`, a `ctJuliaFit` with the draws in
`estimate$rawposterior` and the chain diagnostics in `$sample`: split
R-hat and effective sample size per parameter, divergences, tree depths,
step sizes, E-BFMI (NUTS only), `sampler`, `converged`, `diagnosis`, and
`target`, which names the posterior sampled. This is what
`ctFit(optimize = FALSE)` returns too.

## Details

Every method except `'sample'` writes approximate draws to
`$rawposterior`, drawn from a covariance fitted around the optimum, and
keeps the point estimate at the optimum. The exception is a method run
against a more exact likelihood (`intoverpop = 'quadrature'` or
`filter = 'particle'`), which moves the estimate to the corrected point.
`method = 'sample'` draws from the posterior itself by MCMC and sets
`$estimate$raw` to the posterior mean, as `ctFit(optimize = FALSE)`
does. Both kinds of draws are stored in the same place. A fit with an
approximate uncertainty carries `$uncertainty$settings`, and a sampled
fit carries `$sample`. Methods other than `'sample'` refuse a fit that
has been sampled. `'sample'` can be run again on a sampled fit to draw
more.

To change only the number of draws, use `method = 'stored'`, which
redraws from the covariance the fit already carries.

## Functions

- `ctOptimUncertainty()`: The name of this function in ctsem 3.11.1,
  kept so that older code runs as before. It takes the argument names
  `uncertainty` and `draws`, and has no `intoverpop`, `filter`, `checks`
  or `'sample'`. New code should use `ctFitUncertainty`.

## Methods

- `'checked'`:

  Starts from the normal approximation given by the Hessian, checks it
  against the gradient of the log posterior at draws from that normal,
  and importance samples the directions in which the posterior departs
  from it. A direction is sampled when it would make some parameter's
  95% interval cover less than 90% or more than 99%, or put either tail
  outside 0.5 to 5%. The normal is centred at the posterior mean the
  gradients point to, which differs from the estimate when the posterior
  is skewed. The estimate itself stays as fitted, and
  [`summary()`](https://rdrr.io/r/base/summary.html) reports it beside
  the draws. If nothing needs sampling, the checked widths replace the
  Hessian's where they differ by 5% or more. The sampled draws are
  checked against the gradient (Stein's identity), and a warning names
  any parameter whose mean still looks off by half a posterior sd or
  more, or whose spread looks under half or over twice the posterior's.

  If the posterior follows a curved ridge that sampling cannot correct,
  such as a loading vector whose length is well determined and whose
  direction is not, the Hessian's normal is kept, with a warning naming
  the parameters affected. In models with random effects, the intervals
  of a population sd and of the effects it scales can stay too narrow
  without any warning, so use `'sample'` for intervals to report. The
  check is recorded in `$uncertainty$details$normalcheck`, and any
  sampling in `$uncertainty$details$importance_sampling`.

- `'hessian'`:

  The default, as in
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md): a
  normal approximation from the Hessian at the estimate. The julia
  backend computes the Hessian exactly, and stan by finite differences
  of the gradient. On julia the interval check (`checks`) is made
  alongside and reported, and the normal is left as it is.

- `'bootstrap'`, `'sandwich'`, `'opg'`:

  Built from each subject's contribution to the score: `'bootstrap'` is
  a one-step score bootstrap using the Hessian, `'sandwich'` the
  sandwich estimator (Hessian bread, score covariance meat), and `'opg'`
  the outer product of the scores. With a single subject the
  contributions of single observations are used, with a warning. Fewer
  than ten subjects, or no more score rows than parameters, also give a
  warning.

- `'fullbootstrap'`:

  Resamples subjects and re-optimises each resample, starting from the
  estimate. Stan only.

- `'is'`:

  Adaptive importance sampling against the fit's own objective. Each raw
  parameter is first profiled on both sides of the estimate, from 3 to
  96 standard errors, and a proposal component is added wherever the
  profile shows a tail the curvature misses, such as the skewed tail of
  a variance estimated from few subjects. Rounds of `isitersize` draws
  from these multivariate t components then add components where the
  weights are largest, until the effective sample size reaches `isESS`
  and the Pareto k of the weights is below 0.7 (k needs the loo
  package). The effective sample size, k, the number of evaluations and
  the parameters with a tail are recorded in
  `fit$uncertainty$details$importance_sampling`. A warning is given if k
  is above 0.7, or if a parameter's profile has not dropped off by 96
  standard errors, which suggests an improper posterior. If the Hessian
  covariance is rank deficient, only the directions with curvature are
  sampled, and the others are held at the estimate.

  On a fit whose random effects are integrated (`intoverpop = 'laplace'`
  or `'augmented'`), `'is'` samples the approximate marginal posterior
  the fit used. It corrects the shape of the normal approximation, but
  not any error in the integration itself, such as the Laplace
  approximation's for random effects on variances, in nonlinear dynamics
  or with non-Gaussian indicators. It usually takes several times as
  long as the fit.

  With `intoverpop = 'quadrature'` (julia, Laplace and variational
  fits), `'is'` and `'checked'` draw from a more exact density, in which
  the subjects or groups flagged by the random-effects check are
  integrated by quadrature (`control$nodes` per effect, 5 by default)
  and the rest by the Laplace approximation. The quadrature is made fine
  enough that its error moves no estimate by more than about
  `control$quadratureTolerance` standard errors (0.5 by default), with
  each subject's or group's error counted with its sign. A smaller value
  gives more accurate estimates and takes longer: with thousands of
  subjects, 0.1 took several times as long. The error of this density,
  and how far it can move the estimates, is reported in
  `fit$laplace$correction$error` and in
  [`summary()`](https://rdrr.io/r/base/summary.html). A fit that was not
  refined when it was fitted is first re-optimised on this density,
  unless `correctEstimate = FALSE`. With `'quadrature'`, `'hessian'`
  gives the Hessian of the refined fit at its optimum, as
  `ctFit(intoverpop = 'quadrature')` does. If the refinement stops short
  of its optimum, the Laplace fit's intervals are used, moved to the new
  estimate.

  With `filter = 'particle'`, on fits whose random effects are not
  integrated by Laplace, `'is'` replaces the Kalman filter's likelihood
  by a particle filter's, as computed by
  [`ctParticleLik`](https://cdriveraus.github.io/ctsem/reference/ctParticleLik.md).
  Its settings are entries of `control`: `particles` (1000), `substeps`
  (10), `transition` (`'exponential'` or `'euler'`) and `seed` (1, with
  draw `i` using `seed + i - 1`). The proposal starts from the fit's
  covariance widened by `imisScaleInit` (1.5) and runs at most
  `imisMaxIter` (10) rounds of `isitersize` draws. The estimate moves to
  the weighted posterior mean unless `correctEstimate = FALSE`. The
  result is recorded in `fit$laplace$correction`, or for a particle
  filter in `fit$particle_correction`. If few effective draws remain,
  the approximation is too far from the target to correct by
  reweighting, and `method = 'sample'` with `intoverpop = 'none'` should
  be used instead.

- `'stored'`:

  Computes nothing new. Redraws `finishsamples` normal draws from the
  covariance already on the fit, leaves the covariance, the Hessian and
  the recorded `method` unchanged, and sets
  `$uncertainty$settings$redrawn = TRUE`. Gives a warning if the
  existing draws came from `'is'` or a bootstrap, which normal draws
  from the covariance do not reproduce.

- `'sample'`:

  Markov chain Monte Carlo, described in the next section. Julia backend
  only.

Methods that produce draws aim for an effective sample size of 200
unless told otherwise: `minESS` for `'sample'`, which applies to every
parameter, and `isESS` for `'is'`. Each gives a warning if it ends short
of its target.

## Sampling (`method = 'sample'`)

Draws from the posterior by Markov chain Monte Carlo, as
`ctFit(backend = 'julia', optimize = FALSE)` does and with the same
defaults. The settings are entries of `control`, and `finishsamples` is
not used. The main ones are `chains` (4), `warmup` (200, not used by the
map sampler), `draws` (800, after warmup), `seed`, `saveEffects`
(`FALSE`) and `processes` (`FALSE`). `saveEffects = TRUE` keeps every
draw of every random effect in `fit$sample$effects`. Their means and sds
are always kept, in `fit$sample$effect_mean` and `fit$sample$effect_sd`.
Random effects are reported as each subject's or group's deviation from
the population value, on the raw scale.

*What is sampled* follows the fit's `intoverpop`, as for
`ctFit(optimize = FALSE)`. For `'augmented'` it is the marginal
posterior of the population parameters with the random effects
integrated by the Kalman filter, and for a fit that asked for
`'laplace'` the marginal posterior with the random effects integrated by
the Laplace approximation. For `'none'`, or where `intoverpop = 'auto'`
chose the Laplace approximation, it is the joint posterior of the
population parameters and every subject's random effects. The
`intoverpop` argument changes this: `'none'` for the joint posterior, or
`'laplace'` for the Laplace marginal of a Laplace fit. `'augmented'`
needs a fit made that way, and an augmented fit can only be sampled that
way. Missing values of time independent predictors are sampled too,
unless `tipredMissing` is `'impute'` in the fit's `optimcontrol` or in
`control`.

*The sampler.* `control$sampler` chooses it. For the joint posterior the
default is `'map'`, NUTS through a transport map learnt during warmup,
and `'nuts'` runs NUTS on the parameters and random effects together,
which usually mixes poorly. Both integrate out exactly the random
effects the Kalman filter can (on intercepts, initial means, manifest
means and time dependent predictor effects, in a linear Gaussian model),
and report them drawn from their exact conditional distribution.
`'saem'` draws every random effect by Metropolis-within-Gibbs and the
parameters given them by NUTS. All three sample the same posterior, and
`'map'` is usually the fastest to reach the effective sample size
target.

A marginal posterior is sampled by tempered sequential Monte Carlo
(`'smc'`) where the interval check found that the posterior departs from
the Hessian's normal approximation, and otherwise by NUTS (`'nuts'`),
whose metric is set from the checked widths. If the fit has no interval
check, one is made first. Setting `chains`, `warmup`, `stepsize`,
`adapt_metric`, `adapt_effects` or `settleTol` selects NUTS, and
`sampler = 'smc'` together with `chains`, `warmup` or `stepsize` is an
error. Tempered SMC moves `particles` (200) draws from the normal
approximation to the posterior in stages, then runs as many further
chains as the effective sample size target needs. It is for posteriors
that NUTS, started from the normal approximation, does not converge on,
and costs about the same as NUTS. It reports no R-hat, and stops on the
effective sample size alone. On a marginal posterior, both SMC and NUTS
evaluate a weighted subset of the subjects for most steps and correct
with the full density, so the draws are exact. `coreset = FALSE`
evaluates every subject throughout.

`adapt_metric`, `adapt_effects` and `settleTol` apply to NUTS's warmup.
All three are refused with `'saem'` and `'smc'`, and the last two with
`'map'`. `control$placement` sets where the chains start. `'fit'` starts
them around the fit's estimate, and is the default. If the
random-effects check finds that the Laplace approximation moves an
estimate by a standard error or more, the joint posterior can lie far
from that estimate, so the fit is first refined as with
`intoverpop = 'quadrature'` and the chains start at the refined
estimate. A fit with no check on record is checked first. Leaving
`'randomeffects'` out of `checks` starts the chains at the Laplace
estimate. `'saem'` is for `sampler = 'saem'`, and its default there: it
runs SAEM from the fit's estimate and starts each chain from its
estimate and its draws of the random effects.

*When it stops.* A run stops once every parameter's effective sample
size reaches `minESS` (200) and R-hat is below `rhatTarget` (1.01),
within a budget of four times `draws` per chain, or `maxDraws` if given.
`meanESS` adds a target for the average effective sample size, and
`minESS = 0` takes exactly `draws`. A run that uses its budget without
reaching the targets gives a warning. `settleTol` ends warmup early once
the metric stops changing (off by default). Lower targets give a
quicker, rougher posterior. As a guide, the Monte Carlo error of a
reported quantile is about `1.96 * 2.1 / sqrt(ESS)` posterior sds at the
5% and 95% points and `1.96 * 1.25 / sqrt(ESS)` at the median, and
chains whose means differ by `d` posterior sds give an R-hat of about
`sqrt(1 + d^2)`. With `minESS = 50` and `rhatTarget = 1.07` the error is
about half a posterior sd at the median and more in the tails, and with
the defaults about a fifth of a posterior sd at the median and two
fifths in the tails. In short runs, R-hat is usually the target that
decides when the run stops.

*Processes.* By default the chains run together in this session and
share its `cores`. With `processes = TRUE` they run in up to `cores`
worker processes, each starting its own Julia session, which takes
longer than a short run saves. Each worker runs its share of the chains
in turn, so four chains on two cores run two at a time. The worker
processes stay open for the next call, and
[`ctJuliaKill`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md)`('workers')`
ends them.

## Backend differences

The methods are the same on both backends, and the covariances they give
are comparable, but some details differ.

*Where the result is stored.* A `ctStanFit` gets `fit$stanfit$cov`,
`fit$stanfit$rawposterior` and `fit$stanfit$uncertainty`, and a
`ctJuliaFit` gets `fit$estimate$cov`, `fit$estimate$se`,
`fit$estimate$rawposterior` and `fit$uncertainty`. Both record the
settings used in `$uncertainty$settings`. The draws, covariance and
standard errors are labelled by raw parameter name, as is `coef(fit)`,
so `fit$estimate$se['drift']` or
`apply(fit$estimate$rawposterior, 2, quantile)` can be read by name. The
population sd and correlation parameters are named differently on the
two backends (`popsd_x` and `rawcor_x__y` on julia, `x_SD` and
`y_x_corr` on stan), and come in the same order.

*The Hessian.* Stan takes finite differences of its gradient with one
step size (`control$hessianStep`). Julia computes the Hessian exactly,
and `control$analyticHessian = FALSE` uses finite differences instead,
for comparison. The two agree to about 1e-4 on a well conditioned model.

*Directions with no curvature.* A direction whose curvature is
negligible compared with the largest is left out when the Hessian is
inverted, and so is a direction along which the likelihood is flat
(changing by less than 1.92 over four raw units). Parameters along such
a direction get no standard error and are named in
`fit$identifiability`, and `$uncertainty$cov` records which directions
were left out.

*Whether the intervals are as wide as the curvature allows.* A julia fit
carries `fit$uncertainty$intervalcheck`, which compares each raw
parameter's standard error with `1/sqrt(information[i,i])`, the width
its own curvature gives. The ratio is 1 for a parameter separate from
the rest, and grows as the parameter becomes entangled with others.
Above about 100, the reported width reflects that entanglement more than
the data, and will differ between runs. `$nullmass` gives each
parameter's share of the directions left out of the inversion, and
`$unidentified` names the parameters with enough of it to have an
infinite variance. [`summary()`](https://rdrr.io/r/base/summary.html)
reports their sd, interval and z as `NA`.

*Backend specific arguments.* `method = 'fullbootstrap'`,
`control$bootstrapFitCores`, `control$bootstrapTol` and
`control$parsteps` are for stan only, and are refused on a `ctJuliaFit`.
`method = 'sample'` and `control$analyticHessian` are for julia only.
`cores` sets the number of R worker processes on stan and of engine
threads on julia. The default importance sampling proposals also differ:
`imisScaleInit` and `imisTailScale` are 1.1 on stan and for `'checked'`
on julia, and 1.5 and 1.2 for `'is'` on julia. Set them explicitly to
compare the backends. Summaries of transformed parameters are updated on
stan but not on julia. The raw parameter covariance and draws are
complete on both.

## See also

[`ctFitAddSamples`](https://cdriveraus.github.io/ctsem/reference/ctFitAddSamples.md),
the deprecated stan only predecessor of `method = 'stored'`.
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with
`optimize = FALSE` samples when fitting.
[`ctLaplaceCheck`](https://cdriveraus.github.io/ctsem/reference/ctLaplaceCheck.md)
and
[`ctParticleLik`](https://cdriveraus.github.io/ctsem/reference/ctParticleLik.md)
measure the approximations a julia fit makes, which
`intoverpop = 'quadrature'` and `filter = 'particle'` correct.

## Examples

``` r
# Redraw from the covariance already on the fit, evaluating no model:
fit <- ctFitUncertainty(ctstantestfit, method = 'stored',
  finishsamples = 200)
fit$stanfit$uncertainty$settings$redrawn

# \donttest{
# The Hessian's normal, checked, and importance sampled where it is wrong
fit <- ctFitUncertainty(ctstantestfit, method = 'checked')
summary(fit)$popmeans

# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1))
  model$pars$indvarying <- FALSE
  jfit <- ctFit(ctsemtestdat, model, backend = 'julia', cores = 1,
    priors = TRUE)
  # Posterior draws by MCMC, with lower targets for a quick look
  jfit <- ctFitUncertainty(jfit, method = 'sample', cores = 1,
    control = list(chains = 2, minESS = 50, rhatTarget = 1.05))
  print(jfit$sample)
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
