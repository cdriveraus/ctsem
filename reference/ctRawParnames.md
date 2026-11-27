# ctRawParnames

Gets internal stan parameter names of a ctStanFit object sampled via
stan based on specified substrings. `ctStanParnames` is maintained as a
backward-compatible alias.

## Usage

``` r
ctRawParnames(x, substrings = c("pop_", "popsd"))
```

## Arguments

- x:

  ctStanFit object, sampled
  (`ctFit(backend = 'stan', optimize = FALSE)`).

- substrings:

  vector of character strings, parameter names of the stan model
  containing any of these strings will be returned. Useful strings may
  be 'pop\_' for population means, 'popsd' for population standard
  deviations, or specific combinations such as 'pop_DRIFT' for the
  population means of temporal dynamics parameters

## Value

vector of character strings.

## Examples

``` r
# \donttest{
sunspots<-sunspot.year
sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
id <- 1
time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
ssmodel <- ctModel(type='ct', n.latent=2, n.manifest=1,
 manifestNames='sunspots',
 latentNames=c('ss_level', 'ss_velocity'),
 LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
 DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
 MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
 MANIFESTVAR=diag(0,1), #As per original spec
 CINT=matrix(c(0, 0), nrow=2, ncol=1),
 DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))

#fit
ssfit <- ctFit(datalong, ssmodel, optimize=FALSE, priors=TRUE,
  sampleControl=list(iter=2, chains=1), backend='stan')
ctRawParnames(ssfit,substrings=c('pop_','popsd'))
# }
```
