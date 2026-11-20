# Converts absolute times to intervals for wide format ctsem panel data

Converts absolute times to intervals for wide format ctsem panel data

## Usage

``` r
ctIntervalise(
  datawide,
  Tpoints,
  n.manifest,
  n.TDpred = 0,
  n.TIpred = 0,
  imputedefs = F,
  manifestNames = "auto",
  TDpredNames = "auto",
  TIpredNames = "auto",
  digits = 5,
  mininterval = 0.001,
  individualRelativeTime = TRUE,
  startoffset = 0
)
```

## Arguments

- datawide:

  Wide format data, containing absolute time measurements, to convert to
  interval time scale. See
  [`ctLongToWide`](https://cdriveraus.github.io/ctsem/reference/ctLongToWide.md)
  to easily convert long format data.

- Tpoints:

  Maximum number of discrete time points (waves of data, or measurement
  occasions) for an individual in the input data structure.

- n.manifest:

  number of manifest variables per time point in the data.

- n.TDpred:

  number of time dependent predictors in the data structure.

- n.TIpred:

  number of time independent predictors in the data structure.

- imputedefs:

  if TRUE, impute time intervals based on the measurement occasion (i.e.
  column) they are in, if FALSE (default), set related observations to
  NA. FALSE is recommended unless you are certain that the imputed value
  (mean of the relevant time column) is appropriate. Noise and bias in
  estimates will result if wrongly set to TRUE.

- manifestNames:

  vector of character strings giving variable names of manifest
  indicator variables (without \_Tx suffix for measurement occasion).

- TDpredNames:

  vector of character strings giving variable names of time dependent
  predictor variables (without \_Tx suffix for measurement occasion).

- TIpredNames:

  vector of character strings giving variable names of time independent
  predictor variables.

- digits:

  How many digits to round to for interval calculations.

- mininterval:

  set to lower than any possible observed measurement interval, but
  above 0 - this is used for filling NA values where necessary and has
  no impact on estimates when set in the correct range. (If all observed
  intervals are greater than 1, mininterval=1 may be a good choice)

- individualRelativeTime:

  if TRUE (default), the first measurement for each individual is
  assumed to be taken at time 0, and all other times are adjusted
  accordingly. If FALSE, new columns for an initial wave are created,
  consisting only of observations which occurred at the earliest
  observation time of the entire sample.

- startoffset:

  if 0 (default) uses earliest observation as start time. If greater
  than 0, all first observations are NA, with distance of startoffset to
  first recorded observation.

## Value

The wide format data as a matrix, with the absolute time columns
replaced by the intervals preceding each time point, `dT1`, `dT2`, ...,
and columns named as
[`ctWideNames`](https://cdriveraus.github.io/ctsem/reference/ctWideNames.md)
names them.

## Details

Time column must be numeric! Columns are read in ctsem's wide order, as
[`ctLongToWide`](https://cdriveraus.github.io/ctsem/reference/ctLongToWide.md)
returns them: the manifest variables and then the time dependent
predictors, each by occasion as
[`ctWideNames`](https://cdriveraus.github.io/ctsem/reference/ctWideNames.md)
lists them, then the absolute times `T0`, `T1`, ..., then the time
independent predictors. Columns with exactly those names in another
order are put in that order by name, with a message. Columns with other
names are read by position. An occasion with no time has its
observations set to `NA` unless `imputedefs = TRUE`.

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
