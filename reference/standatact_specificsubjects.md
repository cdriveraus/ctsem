# Adjust standata from ctsem to only use specific subjects

Adjust standata from ctsem to only use specific subjects

## Usage

``` r
standatact_specificsubjects(standata, subjects, timestep = NA)
```

## Arguments

- standata:

  standata

- subjects:

  vector of subjects

- timestep:

  ignored at present

## Value

list of updated structure

## Examples

``` r
d <- standatact_specificsubjects(ctstantestfit$standata, 1:2)
```
