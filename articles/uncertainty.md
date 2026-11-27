# Uncertainty in a Fitted ctsem Model

``` r

library(ctsem)
library(ggplot2)
set.seed(12)
```

A fitted ctsem model reports estimates and their uncertainty. How the
uncertainty is computed depends on three settings:

- `optimize`: sample the posterior (or likelihood) by MCMC, or find its
  mode and approximate the posterior around it.
- `priors`: include ctsem’s weakly informative priors on the raw
  parameters, or leave them out.
- `uncertainty`: for an optimized fit, which approximation to use around
  the mode.

This vignette fits one small model in each of these ways, with both
backends (`backend = 'stan'` and `backend = 'julia'`), and compares the
intervals.

## Sampling and optimizing

**Sampling** (`optimize = FALSE`) draws from the posterior by MCMC, with
Stan’s sampler for `backend = 'stan'` and the Julia engine’s own
samplers for `backend = 'julia'`. The uncertainty is the spread of the
draws, and no normal approximation is involved. This is the most
reliable approach and the slowest.

**Optimizing** (`optimize = TRUE`, the default) finds the mode and
approximates the uncertainty from the curvature there. This takes a
single fit, but the approximation is a normal distribution, which
describes a skewed, bounded or multimodal posterior poorly. Variance
parameters estimated from few subjects often have such posteriors.

For a given `priors` setting both approaches target the same density.
The sections below compare them, show the options between a single
Hessian and full sampling, and then show what changes without priors.

## A small model, fit several ways

The model is a single latent process measured by one indicator with a
fixed intercept, and estimates a continuous time intercept, drift and
diffusion. Eight subjects with eight time points each keep the fits
fast, and the sample is small enough that the uncertainty methods give
visibly different answers.

``` r

n_subjects <- 8
tpoints <- 8

gen_model <- ctModel(
  type = "omx",
  Tpoints = tpoints,
  n.latent = 1,
  n.manifest = 1,
  LAMBDA = matrix(1),
  T0MEANS = matrix(0),
  T0VAR = matrix(0.6),
  DRIFT = matrix(-0.35),
  DIFFUSION = matrix(0.5),
  CINT = matrix(0.15),
  MANIFESTMEANS = matrix(0),
  MANIFESTVAR = matrix(0.25)
)

dat <- data.frame(ctGenerate(gen_model, n = n_subjects, dtmean = 0.4,
  logdtsd = 0.1))

fit_model <- ctModel(
  type = "ct",
  Tpoints = tpoints,
  n.latent = 1,
  n.manifest = 1,
  LAMBDA = matrix(1),
  MANIFESTMEANS = matrix(0),
  CINT = matrix("cint"),
  DRIFT = matrix("drift"),
  DIFFUSION = matrix("diffusion")
)
# No individual differences. By default ctModel() lets T0MEANS and CINT vary
# over subjects, which would add population parameters that eight subjects
# can barely estimate.
fit_model$pars$indvarying <- FALSE
```

The two backends store the estimate and the draws in different places
(`fit$stanfit$rawest` and `fit$stanfit$rawposterior` for stan,
`fit$estimate$raw` and `fit$estimate$rawposterior` for julia).
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`ctExtract()`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md)
work with either and report parameters on their natural scale, so the
comparisons below use
[`ctExtract()`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md).

``` r

drift_interval <- function(fit) {
  draws <- ctExtract(fit)$pop_DRIFT[, 1, 1]
  c(mean = mean(draws), quantile(draws, c(.025, .975)))
}
```

### Sampling

``` r

stan_sample <- ctFit(
  datalong = dat, model = fit_model, backend = "stan",
  cores = 1, priors = TRUE, optimize = FALSE,
  sampleControl = list(chains = 1, iter = 400), verbose = 0
)
#> 
#> SAMPLING FOR MODEL 'ctsm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000866 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 8.66 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 400 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 400 [  5%]  (Warmup)
#> Chain 1: Iteration:  40 / 400 [ 10%]  (Warmup)
#> Chain 1: Iteration:  60 / 400 [ 15%]  (Warmup)
#> Chain 1: Iteration:  80 / 400 [ 20%]  (Warmup)
#> Chain 1: Iteration: 100 / 400 [ 25%]  (Warmup)
#> Chain 1: Iteration: 120 / 400 [ 30%]  (Warmup)
#> Chain 1: Iteration: 140 / 400 [ 35%]  (Warmup)
#> Chain 1: Iteration: 160 / 400 [ 40%]  (Warmup)
#> Chain 1: Iteration: 180 / 400 [ 45%]  (Warmup)
#> Chain 1: Iteration: 200 / 400 [ 50%]  (Warmup)
#> Chain 1: Iteration: 201 / 400 [ 50%]  (Sampling)
#> Chain 1: Iteration: 220 / 400 [ 55%]  (Sampling)
#> Chain 1: Iteration: 240 / 400 [ 60%]  (Sampling)
#> Chain 1: Iteration: 260 / 400 [ 65%]  (Sampling)
#> Chain 1: Iteration: 280 / 400 [ 70%]  (Sampling)
#> Chain 1: Iteration: 300 / 400 [ 75%]  (Sampling)
#> Chain 1: Iteration: 320 / 400 [ 80%]  (Sampling)
#> Chain 1: Iteration: 340 / 400 [ 85%]  (Sampling)
#> Chain 1: Iteration: 360 / 400 [ 90%]  (Sampling)
#> Chain 1: Iteration: 380 / 400 [ 95%]  (Sampling)
#> Chain 1: Iteration: 400 / 400 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 2.159 seconds (Warm-up)
#> Chain 1:                1.707 seconds (Sampling)
#> Chain 1:                3.866 seconds (Total)
#> Chain 1:
```

``` r

julia_sample <- ctFit(
  datalong = dat, model = fit_model, backend = "julia",
  cores = 1, priors = TRUE, optimize = FALSE,
  sampleControl = list(chains = 1, iter = 250), verbose = 0
)
```

Each fit uses one short chain, which is enough for a demonstration. For
real work, use more `chains` and `iter` in `sampleControl`. The drift
interval from each:

``` r

rbind(stan = drift_interval(stan_sample), julia = drift_interval(julia_sample))
#>             mean      2.5%       97.5%
#> stan  -0.7430576 -2.139544 -0.08117748
#> julia -0.7274366 -1.917124 -0.03713031
```

### Optimizing

``` r

stan_map <- ctFit(
  datalong = dat, model = fit_model, backend = "stan",
  cores = 1, priors = TRUE, optimize = TRUE,
  optimcontrol = list(stochastic = FALSE, carefulfit = FALSE, finishsamples = 200)
)
```

``` r

julia_map <- ctFit(
  datalong = dat, model = fit_model, backend = "julia",
  cores = 1, priors = TRUE, optimize = TRUE,
  optimcontrol = list(finishsamples = 200), verbose = 0
)
```

``` r

rbind(stan = drift_interval(stan_map), julia = drift_interval(julia_map))
#>             mean      2.5%      97.5%
#> stan  -0.7051033 -1.945421 -0.1221965
#> julia -0.6331217 -1.925522 -0.1271398
```

By default an optimized fit draws its `rawposterior` from a normal
approximation based on the Hessian (`uncertainty = "hessian"`). The next
section shows the alternatives.

## Alternatives to the Hessian

For an optimized fit, `uncertainty` in `optimcontrol` (when fitting) or
`method` in
[`ctFitUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
(after fitting, from the same estimate) chooses how the uncertainty is
computed.
[`ctFitUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
was called
[`ctOptimUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
in earlier versions, and the old name still works.

- `hessian`: the curvature at the optimum, by finite differences on stan
  and exactly on julia, with draws from the matching normal
  distribution. This is the default and the cheapest. On julia the fit
  also checks the normal approximation, as `checked` does, and says if
  the posterior departs from it. It suits model building, and sampling
  is better for intervals you report.
- `checked`: the `hessian` normal, checked against the gradient of the
  log posterior at 16 to 256 draws, and corrected by importance sampling
  in the directions where it is wrong, such as a skewed parameter or a
  variance that changes the width of other parameters. If the correction
  would cost more than the fit, for instance along a curved ridge, the
  `hessian` normal is kept and a warning names the parameters. When
  nothing is wrong the result is the `hessian` result. The corrected
  directions are sampled with the other parameters held near the
  estimate, so intervals can stay too narrow where a population standard
  deviation widens the effects it scales.
- `opg`: the outer product of per-subject scores, with no Hessian. It is
  mainly a diagnostic.
- `sandwich`: the robust covariance, with the Hessian as bread and the
  covariance of per-subject scores as meat. Useful as a check for
  misspecification.
- `bootstrap`: a one-step bootstrap of the per-subject scores, using the
  Hessian. Much cheaper than refitting, but it needs enough subjects to
  be stable. With few subjects, expect a warning and a wider, noisier
  interval.
- `fullbootstrap`: resample subjects and refit each sample. The most
  expensive method, and only available with `backend = 'stan'`.
- `is`: importance sampling, with the Hessian’s normal as the proposal.
  It reweights the draws towards the objective the fit maximised. Here
  that is the exact posterior. With random effects it includes the
  approximation used to integrate them (see below).

``` r

methods <- c("hessian", "checked", "opg", "sandwich", "bootstrap", "is")

gallery <- function(base_fit, backend_label) {
  do.call(rbind, lapply(methods, function(m) {
    updated <- ctFitUncertainty(base_fit, method = m, finishsamples = 200, cores = 1)
    data.frame(backend = backend_label, method = m, t(drift_interval(updated)))
  }))
}

uncertainty_table <- rbind(gallery(stan_map, "stan"), gallery(julia_map, "julia"))
uncertainty_table
#>    backend    method       mean     X2.5.      X97.5.
#> 1     stan   hessian -0.7137669 -1.939412 -0.15688843
#> 2     stan   checked -0.7088383 -2.148909 -0.11671673
#> 3     stan       opg -0.7660505 -2.748983 -0.12723146
#> 4     stan  sandwich -0.6414727 -1.697333 -0.19003571
#> 5     stan bootstrap -0.6341871 -1.663304 -0.19032650
#> 6     stan        is -0.7223915 -1.842045 -0.07990761
#> 7    julia   hessian -0.6197596 -1.696917 -0.09750200
#> 8    julia   checked -0.7231487 -1.754424 -0.14659718
#> 9    julia       opg -0.7395972 -2.677605 -0.06964777
#> 10   julia  sandwich -0.6140917 -1.753718 -0.12569309
#> 11   julia bootstrap -0.6527327 -1.948373 -0.20200530
#> 12   julia        is -0.7342185 -2.042166 -0.06094387
```

``` r

ggplot(uncertainty_table, aes(x = method, y = mean, ymin = X2.5., ymax = X97.5.,
                               colour = backend)) +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  theme_bw() +
  labs(x = NULL, y = "drift, 95% interval")
```

![](data:image/png;base64,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)

With eight subjects, the score based methods (`opg`, `sandwich`,
`bootstrap`) warn that there are too few subjects for them to be
reliable, and `opg` gives much wider intervals than the others.
`fullbootstrap` runs only with stan and is slower, so it is shown
separately:

``` r

stan_fullboot <- ctFitUncertainty(
  stan_map, method = "fullbootstrap", finishsamples = 8, cores = 1,
  control = list(bootstrapTol = 1e-3)
)
drift_interval(stan_fullboot)
#>       mean       2.5%      97.5% 
#> -0.7077579 -1.1257267 -0.5771079
```

### Estimating the states

Everything above uses `intoverstates = TRUE`, the default. A Kalman
filter integrates over the latent states, and the fit’s objective is the
marginal likelihood (or posterior) of the parameters. With
`intoverstates = FALSE` the states are estimated along with the
parameters. This can help with non-Gaussian indicators, because
`backend = 'julia'` then uses the exact observation density, where the
filter uses an approximation.

Optimizing such a model would find the joint mode of parameters and
states, which is not a useful estimate. With an innovation for every
observation, the states can absorb almost any change in the parameters,
and variance parameters are biased downwards, as they are when random
effects are maximised as if they were parameters.
[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)
therefore refuses `optimize = TRUE` with `intoverstates = FALSE`, and
`optimcontrol$estonly = TRUE` returns the mode anyway. For such a fit
only `uncertainty = 'hessian'` is available, because the other methods
use the marginal likelihood, which the fit did not maximise. The Hessian
is taken with the states profiled out, and on a model this small it is
nearly flat:

``` r

count_gen <- ctModel(
  type = "ct", n.latent = 1, n.manifest = 1,
  manifestNames = "c1", latentNames = "eta1", manifesttype = 3L,
  LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(0.7),
  MANIFESTVAR = matrix(0), T0VAR = matrix(1), T0MEANS = matrix(0),
  CINT = matrix(0), MANIFESTMEANS = matrix(1), Tpoints = 6
)
count_dat <- data.frame(ctGenerate(count_gen, n = 8, Tpoints = 6, backend = "julia"))

count_fit_model <- ctModel(
  type = "ct", n.latent = 1, n.manifest = 1,
  manifestNames = "c1", latentNames = "eta1", manifesttype = 3L,
  LAMBDA = matrix(1), MANIFESTVAR = matrix(0), T0MEANS = matrix(0), CINT = matrix(0),
  MANIFESTMEANS = matrix("cmean"), Tpoints = 6
)
count_fit_model$pars$indvarying <- FALSE

state_explicit_opt <- ctFit(
  count_dat, count_fit_model, backend = "julia",
  intoverstates = FALSE, optimize = TRUE, priors = TRUE, verbose = 0,
  optimcontrol = list(estonly = TRUE)
)
is.null(state_explicit_opt$estimate$se)
#> [1] TRUE
```

`manifesttype = 3` is a count indicator with a Poisson distribution.
With eight subjects and six time points, the profiled Hessian gives no
standard errors at all. Sampling with `optimize = FALSE` draws the
states along with the parameters and makes no Gaussian assumption about
them:

``` r

state_explicit_samp <- ctFit(
  count_dat, count_fit_model, backend = "julia",
  intoverstates = FALSE, optimize = FALSE, priors = TRUE,
  sampleControl = list(chains = 1, iter = 120), verbose = 0
)
state_explicit_samp$estimate$se
#> drift_eta1  diff_eta1      cmean T0var_eta1 
#> 0.65994134 0.19869194 0.01559009 0.26374012
dim(state_explicit_samp$estimate$states)
#> [1] 48  1
```

The sampled fit reports standard errors from the draws, and
`estimate$states` holds the sampled state trajectories.

## With random effects

The model so far has no individual differences. When parameters vary
over subjects, the random effects can be integrated out or sampled, and
this decides what an interval describes.

- **Integrated out**: an optimized fit, with `intoverpop = TRUE` or the
  default `'auto'`. Random effects that shift a mean in a model with
  Gaussian indicators are integrated exactly, and others with a Laplace
  approximation. Every `uncertainty` method above, `is` included, then
  describes this integrated objective and its approximation. Each fit
  checks the Laplace approximation and says if it affects the estimates,
  and `intoverpop = 'quadrature'` refines the fit. The interval is still
  a normal approximation at a mode, and a population standard deviation
  estimated from few subjects usually has a skewed posterior.
- **Sampled**: `optimize = FALSE`, or `ctFitUncertainty(fit, 'sample')`
  after an optimized fit. The parameters and each subject’s random
  effects are drawn together from their joint posterior.

Sampling can be quick. It stops when targets for the effective sample
size and R-hat are met, and lower targets give a rougher posterior with
a known precision. A reported 5% or 95% quantile has a Monte Carlo error
of about `1.96 * 2.1 / sqrt(ESS)` posterior standard deviations, and the
median `1.96 * 1.25 / sqrt(ESS)`. Chains whose means differ by `d`
standard deviations show an R-hat of about `sqrt(1 + d^2)`.
`minESS = 50` with `rhatTarget = 1.07` aims at about half a posterior
standard deviation at the centre. On the models we compared, four chains
of 100 to 200 draws were already more accurate than the Laplace normal
approximation, at two to three times the cost of the fit. For runs of
less than about a minute, keep the chains in the R session
(`processes = FALSE`), because starting Julia in a separate process for
each chain takes tens of seconds.

Here a drift that varies over 30 subjects is summarised three ways. It
is a random effect on the dynamics, so `'auto'` integrates it with the
Laplace approximation:

``` r

set.seed(3)
# The same process, with each subject's drift generated as -exp(raw), where
# raw is normal with mean log(0.35) and sd 0.5 (the sdscale field).
re_gen <- ctModel(type = "ct", Tpoints = tpoints, n.latent = 1, n.manifest = 1,
  LAMBDA = matrix(1), T0MEANS = matrix(0), T0VAR = matrix(0.6),
  DRIFT = matrix("drift|-exp(param)|TRUE|0.5"), DIFFUSION = matrix(0.5),
  CINT = matrix(0.15), MANIFESTMEANS = matrix(0), MANIFESTVAR = matrix(0.25))
re_dat <- data.frame(ctGenerate(re_gen, n = 30, dtmean = 0.4, logdtsd = 0.1,
  backend = "julia", values = c(drift = -0.35)))
re_model <- fit_model
re_model$pars$indvarying[re_model$pars$param %in% "drift"] <- TRUE
```

``` r

fit_time <- system.time(re_fit <- ctFit(re_dat, re_model, backend = "julia",
  priors = TRUE, cores = 1, verbose = 0))[["elapsed"]]
re_fit$args$resolved$intoverpop
#> [1] "laplace"
```

``` r

quick_time <- system.time(re_quick <- ctFitUncertainty(re_fit, "sample",
  cores = 1, control = list(chains = 4, draws = 100, minESS = 50,
    rhatTarget = 1.07, processes = FALSE)))[["elapsed"]]
full_time <- system.time(re_full <- ctFitUncertainty(re_fit, "sample",
  cores = 1, control = list(processes = FALSE)))[["elapsed"]]
```

``` r

drift_rows <- function(fit, label, secs) {
  s <- summary(fit)
  cols <- c("mean", "2.5%", "97.5%")
  data.frame(route = label, seconds = round(secs),
    quantity = c("population mean", "population sd"),
    rbind(s$popmeans["drift", cols], s$popsd["drift", cols]),
    check.names = FALSE, row.names = NULL)
}
rbind(drift_rows(re_fit, "Laplace + Hessian", fit_time),
  drift_rows(re_quick, "short sample", fit_time + quick_time),
  drift_rows(re_full, "default sample", fit_time + full_time))
#>               route seconds        quantity   mean   2.5%  97.5%
#> 1 Laplace + Hessian       8 population mean -0.382 -0.930 -0.110
#> 2 Laplace + Hessian       8   population sd  0.448     NA     NA
#> 3      short sample      29 population mean -0.378 -0.686 -0.137
#> 4      short sample      29   population sd  0.193  0.002  0.811
#> 5    default sample      37 population mean -0.383 -0.730 -0.138
#> 6    default sample      37   population sd  0.226  0.002  1.038
```

Each subject’s data say little about its own drift here. The Laplace
fit’s interval comes from the curvature at its mode, and for the
population standard deviation it can give no interval at all, in which
case [`summary()`](https://rdrr.io/r/base/summary.html) reports `NA` and
the fit warns why. Both samples show a wide, skewed posterior that
reaches down to near zero. The default sample is the reference, and the
short one is a rougher version of it. Its diagnostics show how rough:

``` r

rbind(short = c(max_rhat = max(re_quick$sample$rhat), min_ess = min(re_quick$sample$ess)),
  default = c(max_rhat = max(re_full$sample$rhat), min_ess = min(re_full$sample$ess)))
#>         max_rhat  min_ess
#> short   1.018325 343.7915
#> default 1.009234 609.3049
```

## With and without priors

`priors = TRUE` adds ctsem’s priors on the raw parameters to the
objective, and `priors = FALSE`, the default for
[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md),
leaves them out. (`'randomCorr'` puts a prior on random effect
correlations only, so for a model with no random effects, like this one,
it is the same as `FALSE`.) Together with `optimize` that gives four
combinations:

|  | `priors = TRUE` | `priors = FALSE` |
|----|:---|:---|
| `optimize = TRUE` | MAP: posterior mode, curvature around it | MLE: likelihood maximum, curvature around it |
| `optimize = FALSE` | posterior draws | likelihood draws (see below) |

The first three are shown above. `stan_map` and `julia_map` are MAP
fits, an MLE fit is the same call with `priors = FALSE`, and
`stan_sample` and `julia_sample` are posterior draws. The fourth,
`priors = FALSE` with `optimize = FALSE`, is easy to ask for by
accident, and both backends warn about it:

``` r

stan_lik_sample <- ctFit(
  datalong = dat, model = fit_model, backend = "stan",
  cores = 1, priors = FALSE, optimize = FALSE,
  sampleControl = list(chains = 1, iter = 200), verbose = 0
)
#> 
#> SAMPLING FOR MODEL 'ctsm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000626 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 6.26 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 1: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 1: Iteration:  60 / 200 [ 30%]  (Warmup)
#> Chain 1: Iteration:  80 / 200 [ 40%]  (Warmup)
#> Chain 1: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 1: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 1: Iteration: 120 / 200 [ 60%]  (Sampling)
#> Chain 1: Iteration: 140 / 200 [ 70%]  (Sampling)
#> Chain 1: Iteration: 160 / 200 [ 80%]  (Sampling)
#> Chain 1: Iteration: 180 / 200 [ 90%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 8.149 seconds (Warm-up)
#> Chain 1:                27.189 seconds (Sampling)
#> Chain 1:                35.338 seconds (Total)
#> Chain 1:
range(stan_lik_sample$stanfit$rawposterior)
#> [1]  -150.8883 68001.9921
```

``` r

julia_lik_sample <- ctFit(
  datalong = dat, model = fit_model, backend = "julia",
  cores = 1, priors = FALSE, optimize = FALSE,
  sampleControl = list(chains = 1, iter = 150), verbose = 0
)
range(julia_lik_sample$estimate$rawposterior)
#> [1] -335.0266 1332.4693
```

The posterior draws earlier stayed within a couple of units of zero on
the raw scale. Without a prior, the raw parameters are limited only by
what the likelihood identifies, and for a model and sample this small
that is not much. At these settings the julia backend often warns that
some sampled parameters “reached the region where their transform is
flat to machine precision”, where the posterior is improper and the
chains wander. Warnings are hidden in this rendering, but that one is
worth reading in an interactive session. With stan the same problem
shows as a slow warm-up, and the draws are not reliable either.

For a well identified model with enough data, sampling the likelihood
with `priors = FALSE` is a reasonable choice. For a small or weakly
identified model the raw parameters can wander, as they do here, so
`priors = TRUE` is the safer choice for sampling, whatever the default
for optimizing.

## Summary

- For a given `priors` setting, sampling and optimizing target the same
  density. Sampling returns draws, and optimizing returns a point
  estimate with an approximation around it.
- For an optimized fit, `uncertainty` chooses the approximation.
  `hessian` is the cheapest and the default, and `checked` corrects it
  where it is wrong. `opg`, `sandwich` and `bootstrap` use per-subject
  scores and need enough subjects. `fullbootstrap` (stan only) refits
  resampled data and is the most expensive, and `is` uses the Hessian’s
  normal as an importance sampling proposal.
- With random effects, the uncertainty of an optimized fit describes the
  integrated objective, including its approximation, whichever method is
  used. `ctFitUncertainty(fit, 'sample')` samples the joint posterior,
  and a lower `minESS` and higher `rhatTarget` make that quicker, with a
  known precision.
- `intoverstates = FALSE` is for sampling, which draws the states along
  with the parameters. Optimizing is refused unless
  `optimcontrol$estonly = TRUE` asks for the joint mode, and then only
  `hessian` is available.
- Priors change the objective: MAP or MLE when optimizing, the posterior
  or the likelihood when sampling. Sampling the likelihood of a small or
  weakly identified model needs care, because the chains can wander.
