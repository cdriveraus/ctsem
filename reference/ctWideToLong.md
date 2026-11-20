# ctWideToLong Convert ctsem wide to long format

ctWideToLong Convert ctsem wide to long format

## Usage

``` r
ctWideToLong(
  datawide,
  Tpoints,
  n.manifest,
  n.TDpred = 0,
  n.TIpred = 0,
  manifestNames = "auto",
  TDpredNames = "auto",
  TIpredNames = "auto"
)
```

## Arguments

- datawide:

  ctsem wide format data

- Tpoints:

  number of measurement occasions in data

- n.manifest:

  number of manifest variables

- n.TDpred:

  number of time dependent predictors

- n.TIpred:

  number of time independent predictors

- manifestNames:

  Character vector of manifest variable names.

- TDpredNames:

  Character vector of time dependent predictor names.

- TIpredNames:

  Character vector of time independent predictor names.

## Value

A long format matrix, one row per subject and time point, with the id,
the time interval `dT` (absolute times are restored by
[`ctDeintervalise`](https://cdriveraus.github.io/ctsem/reference/ctDeintervalise.md)),
the manifest variables and the predictors.

## Details

Names must account for \*all\* the columns in the data - i.e. do not
leave certain variables out just because you do not need them.

## Examples

``` r
#create wide data
wideexample <- ctLongToWide(datalong = ctstantestdat, id = "id", 
time = "time", manifestNames = c("Y1", "Y2"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2","TI3"))

wide <- ctIntervalise(datawide = wideexample, Tpoints = 10, n.manifest = 2, 
n.TDpred = 1, n.TIpred = 3, manifestNames = c("Y1", "Y2"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2","TI3") )

 #Then convert to long format
 longexample <- ctWideToLong(datawide = wideexample, Tpoints=10, 
 n.manifest=2, manifestNames = c("Y1", "Y2"),
 n.TDpred=1, TDpredNames = "TD1", 
 n.TIpred=3, TIpredNames = c("TI1", "TI2","TI3"))

 #Then convert the time intervals to absolute time
 long <- ctDeintervalise(datalong = longexample, id='id', dT='dT')
 head(long,22)
```
