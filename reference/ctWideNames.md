# Column names for wide ctsem data

Gives the default column names for wide ctsem data. Mostly used inside
ctsem.

## Usage

``` r
ctWideNames(
  n.manifest,
  Tpoints,
  n.TDpred = 0,
  n.TIpred = 0,
  manifestNames = "auto",
  TDpredNames = "auto",
  TIpredNames = "auto"
)
```

## Arguments

- n.manifest:

  number of manifest variables per time point in the data.

- Tpoints:

  Maximum number of discrete time points (waves of data, or measurement
  occasions) for an individual in the input data structure.

- n.TDpred:

  number of time dependent predictors in the data structure.

- n.TIpred:

  number of time independent predictors in the data structure.

- manifestNames:

  vector of character strings giving column names of manifest indicator
  variables

- TDpredNames:

  vector of character strings giving column names of time dependent
  predictor variables

- TIpredNames:

  vector of character strings giving column names of time independent
  predictor variables

## Value

A character vector of column names for a wide format dataset, with the
manifest variables and time dependent predictors at each time point,
then the intervals `dT1`, `dT2`, ..., then the time independent
predictors.

## Examples

``` r
# Column names, in order, for wide data with two variables measured at
# three occasions, one time dependent and one time independent predictor.
# dT1 and dT2 are the time intervals before the second and third occasion.
ctWideNames(n.manifest = 2, Tpoints = 3, n.TDpred = 1, n.TIpred = 1)

# With the variables' own names.
ctWideNames(n.manifest = 2, Tpoints = 3,
  manifestNames = c('LeisureTime', 'Happiness'))
```
