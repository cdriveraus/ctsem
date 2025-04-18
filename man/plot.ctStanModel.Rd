% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.ctStanModel.R
\name{plot.ctStanModel}
\alias{plot.ctStanModel}
\title{Prior plotting}
\usage{
\method{plot}{ctStanModel}(
  x,
  rows = "all",
  wait = FALSE,
  nsamples = 1e+06,
  rawpopsd = "marginalise",
  inddifdevs = c(-1, 1),
  inddifsd = 0.1,
  plot = TRUE,
  ...
)
}
\arguments{
\item{x}{ctStanModel object as generated by \code{\link{ctModel}} with type='ct' or 'dt'.}

\item{rows}{vector of integers denoting which rows of ctstanmodel$pars to plot priors for. 
Character string 'all' plots all rows with parameters to be estimated.}

\item{wait}{If true, user is prompted to continue before plotting next graph.}

\item{nsamples}{Numeric. Higher values increase fidelity (smoothness / accuracy) of density plots, at cost of speed.}

\item{rawpopsd}{Either 'marginalise' to sample from the specified (in the ctstanmodel) 
prior distribution for the raw population standard deviation, or a numeric value to use for the raw population standard deviation
for all subject level prior plots - the plots in dotted blue or red.}

\item{inddifdevs}{numeric vector of length 2, setting the means for the individual differences distributions.}

\item{inddifsd}{numeric, setting the standard deviation of the population means used to generate individual
difference distributions.}

\item{plot}{If FALSE, ouputs list of GGplot objects that can be further modified.}

\item{...}{not used.}
}
\description{
Plots priors for free model parameters in a ctStanModel.
}
\details{
Plotted in black is the prior for the population mean. In red and blue are the subject level priors that result
given that the population mean is estimated as 1 std deviation above the mean of the prior, or 1 std deviation below. 
The distributions around these two points are then obtained by marginalising over the prior for the raw population std deviation - 
so the red and blue distributions do not represent any specific subject level prior, but rather characterise the general amount
and shape of possible subject level priors at the specific points of the population mean prior.
}
\examples{
model <- ctModel(type='ct',
manifestNames='sunspots', 
latentNames=c('ss_level', 'ss_velocity'),
LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
# MANIFESTVAR=matrix(0, nrow=1, ncol=1),
CINT=matrix(c(0, 0), nrow=2, ncol=1),
DIFFUSION=matrix(c(
  0, 0,
  0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
  
plot(model,rows=8)
}
