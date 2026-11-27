# Defunct functions in ctsem

These functions were part of ctsem 3.11.1 and have since been removed.
Each now stops with a message naming its replacement.

- `ctCheckFit`: use
  [`ctPostPredPlots`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md),
  whose `'LaggedCovariance'` and `'MeanTrajectory'` panels replace it,
  or
  [`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
  for the lagged covariance check alone.

- `ctFitUpdateModel` and `ctStanUpdModel`: use
  [`ctFitUpdate`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md),
  which takes new data and new arguments.

- `ctModelCoverage_check`: renamed
  [`ctCoverageCheck`](https://cdriveraus.github.io/ctsem/reference/ctCoverageCheck.md),
  with camelCase arguments (`fittingModel`, `fitArgs`, `fitCores`,
  `generateCores`, `plotEvery`).

## Usage

``` r
ctFitUpdateModel(...)

ctStanUpdModel(...)

ctCheckFit(...)

ctModelCoverage_check(...)
```

## Arguments

- ...:

  Ignored.
