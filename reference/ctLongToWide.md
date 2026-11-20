# ctLongToWide Restructures time series / panel data from long format to wide format for ctsem analysis

ctLongToWide Restructures time series / panel data from long format to
wide format for ctsem analysis

## Usage

``` r
ctLongToWide(
  datalong,
  id,
  time,
  manifestNames,
  TDpredNames = NULL,
  TIpredNames = NULL
)
```

## Arguments

- datalong:

  dataset in long format, including subject/id column, observation time
  (or change in observation time, with 0 for first observation) column,
  indicator (manifest / observed) variables, any time dependent
  predictors, and any time independent predictors.

- id:

  character string giving column name of the subject/id column

- time:

  character string giving column name of the time columnn

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

A wide format matrix, one row per subject: each manifest variable and
time dependent predictor at each time point (`Y1_T0`, `Y1_T1`, ...), the
observation times `T0`, `T1`, ..., then the time independent predictors.

## Details

Time column must be numeric

## See also

[`ctIntervalise`](https://cdriveraus.github.io/ctsem/reference/ctIntervalise.md)

## Examples

``` r
wideexample <- ctLongToWide(datalong = ctstantestdat, id = "id", 
time = "time", manifestNames = c("Y1", "Y2"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2","TI3"))

#Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
wide <- ctIntervalise(datawide = wideexample, Tpoints = 10, n.manifest = 2, 
n.TDpred = 1, n.TIpred = 3, manifestNames = c("Y1", "Y2"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2","TI3") )
 
print(wide)
```
