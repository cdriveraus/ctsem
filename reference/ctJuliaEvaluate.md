# Evaluate a julia model's likelihood and gradient at given parameters

Evaluates the objective that the julia backend optimises, the log
likelihood with the random effects integrated out as `intoverpop`
specifies, and its gradient, at a vector of raw (unconstrained)
parameters. Useful for checking a fit, profiling by hand, or passing
ctsem's likelihood to another optimiser or sampler.

## Usage

``` r
ctJuliaEvaluate(
  object,
  pars = NULL,
  gradient = TRUE,
  contributions = FALSE,
  gradient_method = c("adjoint", "forward"),
  cores = getOption("mc.cores", 2L)
)
```

## Arguments

- object:

  A fit from `ctFit(..., backend = 'julia')`, or a model prepared
  without fitting by `ctFit(..., backend = 'julia', fit = FALSE)` (class
  `ctJuliaModel`).

- pars:

  Raw (unconstrained) parameters, in the order of `coef(fit)`. Defaults
  to the fit's estimate, and must be given for a prepared model.

- gradient:

  Also return the gradient with respect to `pars`.

- contributions:

  Return subject and row contribution details when available.

- gradient_method:

  Either "adjoint" (reverse mode, the default) or "forward"
  (ForwardDiff). Both compute the same gradient. The cost of "adjoint"
  does not grow with the number of parameters, so it is faster, and more
  so for larger models. See `optimcontrol$gradient` in
  [`ctOptimControl`](https://cdriveraus.github.io/ctsem/reference/ctOptimControl.md).

- cores:

  The most engine threads this evaluation uses. `NA` keeps the limit
  already set.

## Value

A list containing the log likelihood and, when requested, the gradient
and the contributions.

## Examples

``` r
# \donttest{
# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1))
  model$pars$indvarying <- FALSE
  prepared <- ctFit(ctstantestdat, model, backend = 'julia', fit = FALSE)
  start <- ctJuliaEvaluate(prepared, pars = rep(0, 6), cores = 1)
  str(start)
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
