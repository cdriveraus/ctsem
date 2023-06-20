
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->
[![R-CMD-check](https://github.com/cdriveraus/ctsem/workflows/R-CMD-check/badge.svg)](https://github.com/cdriveraus/ctsem/actions)
<!-- badges: end -->

**See the NEWS file for recent updates, and below for quick start\!**

ctsem allows for easy specification and fitting of a range of continuous
and discrete time dynamic models, including multiple indicators (dynamic
factor analysis), multiple, potentially higher order processes, and time
dependent (varying within subject) and time independent (not varying
within subject) covariates. Classic longitudinal models like latent
growth curves and latent change score models are also possible. Version
1 of ctsem provided SEM based functionality by linking to the OpenMx
software, allowing mixed effects models (random means but fixed
regression and variance parameters) for multiple subjects. For version 2
of the R package ctsem, we include a hierarchical specification and
fitting routine that uses the Stan probabilistic programming language,
via the rstan package in R. This allows for all parameters of the
dynamic model to individually vary, using an estimated population mean
and variance, and any time independent covariate effects, as a prior.
Version 3 allows for state dependencies in the parameter specification
(i.e. time varying parameters).

The current manual is at
<https://cran.r-project.org/package=ctsem/vignettes/hierarchicalmanual.pdf>.
The original ctsem is documented in a JSS publication (Driver, Voelkle,
Oud, 2017), and in R vignette form at
<https://cran.r-project.org/package=ctsemOMX/vignettes/ctsem.pdf>,
however these OpenMx based functions have been split off into a sub
package, ctsemOMX. For most use cases the newer formulation (with Kalman
filtering coded in Stan) is faster, more robust, and more flexible, and
both default to maximum likelihood. For cases with many subjects, few
time points, and no individual differences in timing, ctsemOMX may be
faster.

For questions (or to see past answers) please use
<https://github.com/cdriveraus/ctsem/discussions>

To cite ctsem please use the citation(“ctsem”) command in R.

### To install the github version, first install rstan and Rtools, then from a fresh R session:

``` r
remotes::install_github('cdriveraus/ctsem', INSTALL_opts = "--no-multiarch", dependencies = c("Depends", "Imports"))
```

### Or just use the CRAN version, but rstan compiler setup is needed separately for some models:

``` r
install.packages('ctsem')
```

### Troubleshooting Rstan / Rtools install for Windows:

Ensure recent version of R and Rtools is installed. If the
installctsem.R code has never been run before, be sure to run that (see
above).

Place this line in \~/.R/makevars.win , and if there are other lines,
delete them:

    CXX14FLAGS += -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations

see  for details

In case of compile errors like `g++ not found`, ensure the devtools
package is installed:

``` r
install.packages('devtools')
```

and include the following in your .Rprofile, replacing c:/Rtools with
the appropriate path – sometimes Rbuildtools/4.0/ .

``` r
library(devtools)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(PATH = paste("C:/Rtools/mingw_64/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
```

### Quick start – univariate panel data with covariate effects on parameters

\#’ The basic long data structure. Diet, (our covariate) is a
categorical variable so needs dummy / ‘one hot’ encoding.

``` r
head(ChickWeight) 
```

\#’ Setup dummy coding

``` r
library(data.table)
library(mltools)
chickdata <- one_hot(as.data.table(ChickWeight),cols = 'Diet')
```

\#’ Scaling of continuous variables makes for easier estimation and more
sensible default priors (if used). Time intervals can also benefit

``` r
chickdata$weight <- scale(chickdata$weight) 
head(chickdata) #now we have the four diet categories
```

\#’ Setup continuous time model – in this case we are estimating a
regular first order autoregressive

``` r
library(ctsem)

m <- ctModel(
  LAMBDA=diag(1), #Factor loading matrix of latent processes on measurements, fixed to 1
  type = 'stanct', #Could specify 'standt' here for discrete time.
  tipredDefault = FALSE, #limit covariate effects on parameters to those explicitly specified
  manifestNames='weight', #Observed measurements of the latent processes
  latentNames='Lweight', #Names here simply make parameters and plots more interpretable
  TIpredNames = paste0('Diet_',2:4), #Covariates, in this case one category needs to be baseline...
  DRIFT='a11 | param', #normally self feedback (diagonal drift terms) are restricted to negative
  MANIFESTMEANS=0, #For identification CINT is normally zero with this freely estimated
  CINT='cint ||||Diet_2,Diet_3,Diet_4', #diet covariates specified in 5th 'slot' (four '|' separators)
  time='Time',
  id='Chick')
```

\#’ View model in pdf/ latex form

``` r
ctModelLatex(m)
```

\#’ Fit model to data – here using priors because Hessian problems are
reported otherwise

``` r
f <- ctStanFit(chickdata,m,priors=TRUE) 
```

\#’ Summarise fit, view covariate effects – Diets 3 and 4 seem most
obviously successful

``` r
s=summary(f)

print(s$tipreds )
```

\#’ Predictions conditional on all earlier data

``` r
ctKalman(f,plot=TRUE,subjects=2:4,kalmanvec=c('yprior','ysmooth')) 
```

\#’ Predictions conditional only on covariates, showing 1 chick from
each diet

``` r
ctKalman(f,plot=T, 
  subjects=as.numeric(chickdata$Chick[!duplicated(ChickWeight$Diet)]),
  removeObs = T,polygonalpha=0)
```

\#’ Plot temporal regression coefficients conditional on time interval –
increases in this case\!

``` r
ctStanDiscretePars(f,plot=T) 
```

\#’ Other useful functions:

\#’ Compare two fits: ctChisqTest()

\#’ Fit and summarise / plot a list of models: ctFitMultiModel()

\#’ Add samples to fit to increase estimate precision: ctAddSamples()

\#’ Return dynamic system parameters in matrix forms:
ctStanContinuousPars()

\#’ Compute cross validation statistics: ctLOO()

\#’ Plot time independent predictor (covariate effects on parameters):
ctStanTIpredEffects()

\#’ Generate data from a specified model of fixed parameters:
ctGenerate()

\#’ Generate data from a specified model of fixed and free parameters /
priors: ctStanGenerate()

\#’ Generate data from a fitted model: ctStanGenerateFromFit()

\#’ Get samples from the fitted object: ctExtract()

\#’ In samples, pop\_DRIFT refers to the population drift matrix,
subj\_DRIFT refers to the subject matrix. Subject matrices only computed
for max likelihood / posterior mode by default, and found in the
\(stanfit\)transformedparsfull object.
