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
