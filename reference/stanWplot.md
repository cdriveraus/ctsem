# Runs stan, and plots sampling information while sampling.

Runs stan, and plots sampling information while sampling.

## Usage

``` r
stanWplot(object, iter = 2000, chains = 4, ...)
```

## Arguments

- object:

  stan model object

- iter:

  Number of iterations

- chains:

  Number of chains

- ...:

  All the other regular arguments to stan()

## Value

The stanfit object `rstan::sampling` returns.

## Details

On windows, requires Rtools installed and able to be found by
pkgbuild::rtools_path()

## Examples

``` r
library(rstan)
#### example 1 
scode <- "
parameters {
  real y[2]; 
} 
model {
  y[1] ~ normal(0, .5);
  y[2] ~ double_exponential(0, 2);
} 
"
#Uncomment the following lines -- launches rscript not compatible with cran check.
#sm <- stan_model(model_code = scode)
#fit1 <- stanWplot(object = sm,iter = 100000,chains=2,cores=1)
```
