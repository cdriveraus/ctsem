# Check interval coverage by simulation

Performs a coverage check analysis by generating data from a model,
fitting it multiple times with different fit arguments, and plotting the
results. Works with `fitArgs`/ `fittingModel` combinations that fit with
either `backend='julia'` (the
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)
default) or `backend='stan'`.

## Usage

``` r
ctCoverageCheck(
  initialData,
  fittingModel,
  niter,
  fitArgs,
  cores = getOption("mc.cores", 2L),
  fitCores = 1,
  generateCores = fitCores,
  plotEvery = max(c(10, cores))
)
```

## Arguments

- initialData:

  An initial dataset to fit to determine 'true' parameters for further
  generation.

- fittingModel:

  A ctModel object used for fitting the data

- niter:

  Number of iterations to run

- fitArgs:

  Named list of fit argument sets to test (e.g., list(boot =
  list(optimcontrol = list(uncertainty = 'bootstrap')), hess =
  list(optimcontrol = list(uncertainty = 'hessian'))))

- cores:

  The most CPU cores the check uses at once. Iterations run
  `cores %/% fitCores` at a time, each fit on `fitCores`.

- fitCores:

  Number of cores for each inner
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) call,
  at most `cores`. It replaces any `cores` given in `fitArgs`.

- generateCores:

  Number of cores for
  [`ctGenerateFromFit`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromFit.md)
  when generating replicated datasets, at most `cores`. Defaults to
  `fitCores`; set it explicitly to use more cores for generation.

- plotEvery:

  Print plots every n iterations; by default every 10, or every `cores`
  iterations when that is more.

## Value

A list containing the results data.table and final plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Does the Hessian interval cover as often as it says? 100 datasets
# generated from the fit to ctstantestdat, each fitted twice.
model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
  LAMBDA = matrix(1))
check <- ctCoverageCheck(ctstantestdat, model, niter = 100,
  fitArgs = list(
    hessian = list(optimcontrol = list(uncertainty = 'hessian')),
    bootstrap = list(optimcontrol = list(uncertainty = 'bootstrap'))),
  cores = 4)
} # }
```
