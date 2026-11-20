# Check a julia-backend fit against a particle-filter likelihood

The Kalman filter ctsem fits with approximates the distribution of the
latent states by a Gaussian. It linearises the drift over each
prediction step, and approximates non-Gaussian measurement updates by
Gaussian ones. For a linear model with Gaussian indicators it is exact,
and otherwise its error cannot be seen from the fit alone. This function
computes the same marginal log likelihood at the fit's estimates with a
bootstrap particle filter, which makes none of these approximations, and
reports the difference, in total and per row.

## Usage

``` r
ctParticleLik(
  fit,
  particles = 2000,
  substeps = 20,
  transition = c("exponential", "euler"),
  seed = 1,
  replicates = 1,
  cores = getOption("mc.cores", 2L)
)
```

## Arguments

- fit:

  A fit from `ctFit(..., backend = 'julia')` with the default
  `intoverpop`. Random effects carried as extra states are integrated by
  the particles. Fits using the Laplace approximation are refused.

- particles:

  Particles per subject.

- substeps:

  Transition steps per observation interval, independent of the fit's
  own settings, but never coarser than the fit's `maxtimestep` or
  automatic step size.

- transition:

  `'exponential'` uses the filter's own discretisation at each
  particle's state, which is exact for the linear parts of the model, so
  the comparison shows the error of the filter's Gaussian approximation.
  `'euler'` is plain Euler-Maruyama, which shares no approximation with
  the filter but needs several times more substeps.

- seed:

  Seed for every draw, so the result is reproducible.

- replicates:

  Number of independent runs (with seeds `seed`, `seed + 1`, ...). With
  one run, an approximate standard error is computed from the effective
  sample sizes, and with more, the standard deviation across runs is
  reported.

- cores:

  Engine threads. Each subject has its own random stream, so the result
  does not depend on the number of cores.

## Value

A list: `loglik`, the particle estimate (mean over replicates); `se`;
`fit_loglik`, the filter's log likelihood from the fit; `difference`,
particle minus filter; `rows`, a data frame with the subject, time,
particle and filter log likelihood of every row and their difference;
`replicates`, the individual estimates; `ess_min`, the smallest
effective sample size any row saw; and the settings.

## Details

Each particle moves through the process with drift, diffusion and every
state dependent parameter evaluated at its own state, is weighted by the
density of the row's observations given that state (Gaussian, binary,
ordinal, count and censored indicators are all handled), and is
resampled when the weights become uneven. The estimate is unbiased, and
its Monte Carlo error falls with the square root of `particles`.

It is a check at fixed parameters, and does not estimate them. The
comparison can be trusted when `se` is well below the difference that
matters (raise `particles`, or set `replicates` for an empirical
standard error), and doubling `substeps` does not change the answer.

## See also

[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
with `method = 'is'` and `filter = 'particle'` applies the correction
this measures, to the fit's posterior draws.

## Examples

``` r
# \donttest{
# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1))
  model$pars$indvarying <- FALSE
  fit <- ctFit(ctstantestdat, model, backend = 'julia', cores = 1)
  # A linear Gaussian model, where the filter is exact: the difference
  # should be within a few standard errors of zero.
  check <- ctParticleLik(fit, particles = 1000, substeps = 5)
  check[c('loglik', 'se', 'fit_loglik', 'difference')]
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
