# Discretise long format continuous time (ctsem) data to specific timestep

Rounds timing information so equal intervals, according to specified
timestep, are achieved. Each subject's times are moved to the nearest
point of a grid that starts at the subject's first time. Values that
land on the same point are averaged, with a warning. NA's are inserted
in other columns as necessary. Any columns specified by TDpredNames have
zeroes rather than NA's inserted (because some estimation routines do
not tolerate NA's in covariates). Columns specified by TIpredNames are
left as NA's in inserted rows.

## Usage

``` r
ctDiscretiseData(
  dlong,
  timestep,
  timecol = "time",
  idcol = "id",
  TDpredNames = NULL,
  TIpredNames = NULL
)
```

## Arguments

- dlong:

  Long format data

- timestep:

  Positive real value to discretise

- timecol:

  Name of column containing absolute (not intervals) time information.

- idcol:

  Name of column containing subject id variable.

- TDpredNames:

  Vector of column names of any time dependent predictors

- TIpredNames:

  Vector of column names of any time independent predictors

## Value

long format ctsem data.

## Examples

``` r
long <- ctDiscretiseData(dlong=ctsemtestdat, timestep = .1,
TDpredNames=c('TD1'),TIpredNames=c('TI1','TI2','TI3'))
```
