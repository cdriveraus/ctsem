% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctStanPostPredict.R
\name{ctPostPredData}
\alias{ctPostPredData}
\title{Create a data.table to compare data generated from a ctsem fit with the original data.}
\usage{
ctPostPredData(fit, residuals = F)
}
\arguments{
\item{fit}{A fitted ctsem model.}

\item{residuals}{If set to TRUE, includes residuals in the comparison.}
}
\value{
A data table containing the comparison between generated and original data.
}
\description{
This function allows for easy comparison of data generated from a fitted ctsem model
with the original data used to fit the model. It provides options to include residuals
in the comparison.
}
\examples{
data_comparison <- ctPostPredData(ctstantestfit)

}
\seealso{
Other ctsem functions for model fitting and analysis.
}
