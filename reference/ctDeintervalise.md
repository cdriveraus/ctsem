# ctDeintervalise

Converts intervals in ctsem long format data to absolute time

## Usage

``` r
ctDeintervalise(datalong, id = "id", dT = "dT", startoffset = 0)
```

## Arguments

- datalong:

  data to use, in ctsem long format (attained via function ctWideToLong)

- id:

  character string denoting column of data containing numeric identifier
  for each subject.

- dT:

  character string denoting column of data containing time interval
  preceding observations in that row.

- startoffset:

  Number of units of time to offset by when converting.

## Value

`datalong`, with the interval column `dT` replaced by absolute times and
renamed `time`.

## Examples

``` r
# Wide data, one row per subject, with the time intervals between
# occasions in the columns dT1 to dT5.
data(ctExample1)
long <- ctWideToLong(ctExample1, Tpoints = 6, n.manifest = 2,
  manifestNames = c('LeisureTime', 'Happiness'))

# In long format the dT column holds the interval before each row.
head(long)

# The intervals summed within each subject into a time column, starting
# at 0.
head(ctDeintervalise(long))

# Starting each subject at time 10 instead.
head(ctDeintervalise(long, startoffset = 10))
```
