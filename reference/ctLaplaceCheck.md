# Measure, and optionally correct, the Laplace approximation error in a fit

Recomputes the log marginal likelihood of a `backend='julia'`,
`intoverpop='laplace'` fit by adaptive Gauss-Hermite quadrature over the
random effects, and reports how far the estimate would move if the
random effects were integrated exactly.

## Usage

``` r
ctLaplaceCheck(
  fit,
  nodes = 5L,
  correction = TRUE,
  step = 0.001,
  refine = FALSE,
  maxiter = 50L,
  cores = getOption("mc.cores", 2L),
  verbose = 0L
)
```

## Arguments

- fit:

  A `ctJuliaFit` fitted with `intoverpop='laplace'`.

- nodes:

  Quadrature nodes per random effect. The cost is `nodes^k` likelihood
  evaluations per subject or group, where `k` is the number of random
  effects at that level. 5 is usually enough to locate the maximum, and
  9 to settle the value.

- correction:

  Compute the first-order correction to the estimate. Costs `2 * npar`
  quadrature evaluations and needs the fit's Hessian, so it is skipped
  when the fit has no uncertainty computed.

- step:

  Finite-difference step for the gap gradient, on the raw scale.

- refine:

  Continue the fit on the quadrature objective instead of correcting
  linearly. This runs the refined fit of `intoverpop = 'quadrature'`
  (see [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md))
  from this fit's estimate, and reports the point it reaches in
  `refined`.

- maxiter:

  Iteration limit for `refine`.

- cores:

  Engine threads for the quadrature.

- verbose:

  Integer. 1 or more prints progress.

## Value

A list with `gap` (quadrature minus Laplace log marginal likelihood at
the estimate), `gap_per_subject`, `nodes`, and, when `correction` is
`TRUE`, a `parameters` data frame with one row per raw parameter giving
`estimate`, `delta`, `corrected`, `se` and `delta_se`. With `refine` it
also contains `refined`. `at` says whether the estimate checked was the
Laplace optimum (`'laplace'`) or a refined estimate (`'corrected'`).

`dropped_directions` counts the directions in which the data give little
or no information, those named in `fit$identifiability`. The correction
is not computed along them, and is reported as zero there.

## Details

The quadrature is centred at each subject's random-effect mode and
scaled by the curvature there, both taken from the fit. The curvature's
eigenvalues are floored at one, the curvature of the prior, so that the
quadrature nodes stay within the prior where the likelihood is flat in
some direction.

The gap is the quadrature log likelihood minus the Laplace log
likelihood that the fit maximised. With `nodes = 1` the gap is zero
wherever the curvature is at least the prior's.

The correction is a single Newton step,
`delta = solve(-hessian, gradient of (quadrature - laplace))`, at the
estimate. Compare it with the standard errors. `delta / se` well below
one means the approximation does not limit the estimate. Near or above
one, it does, and the interval may not cover the true value.

To refine the fit with the random effects integrated by quadrature, use
`ctFit(intoverpop = 'quadrature')`, or
`ctFitUncertainty(fit, intoverpop = 'quadrature')` after fitting. On a
refined fit the correction reported here is the further Newton step from
the refined estimate, and `at` in the result is `'corrected'`. A small
further step means the refinement converged.

Nested grouping levels are supported. The quadrature runs over each
level given the levels above it, so the cost grows linearly with the
number of groups.

## See also

[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
with `method = 'sample'` and `intoverpop = 'none'` samples the joint
posterior of parameters and random effects, which avoids the
approximation. It is much slower, and the better choice when `delta_se`
is large. If the fit's random-effects check finds a move of a standard
error or more, it refines the fit first and starts the chains there.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- ctsemtestdat
model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
  LAMBDA = matrix(1))
model$pars$indvarying <- model$pars$matrix %in% 'MANIFESTMEANS'
fit <- ctFit(data, model, backend = 'julia', intoverpop = 'laplace')

check <- ctLaplaceCheck(fit)
check                     # gap, and the largest corrections in standard errors
check$parameters          # the full per-parameter table
} # }
```
