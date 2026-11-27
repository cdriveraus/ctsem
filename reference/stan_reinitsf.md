# Quickly initialise stanfit object from model and data

Quickly initialise stanfit object from model and data

## Usage

``` r
stan_reinitsf(model, data, fast = FALSE)
```

## Arguments

- model:

  stanmodel

- data:

  standata

- fast:

  Use cut down form for speed

## Value

stanfit object

## Examples

``` r
sf <- stan_reinitsf(ctstantestfit$stanmodel,ctstantestfit$standata)
```
