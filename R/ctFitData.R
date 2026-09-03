#' Get the data used to fit a model, in long format, for any backend
#'
#' A stan fit and a julia fit keep the input data differently: a stan fit's
#' \code{$data} is the prepared \code{standata} structure (wide, engine-
#' internal), while a julia fit's \code{$data} is the original long-format
#' data frame. This function returns the same thing -- the original long
#' data frame, one row per observation -- for either backend, so code that
#' works with fitted data does not need to know which backend produced the
#' fit.
#'
#' @param fit A fitted ctsem model, from \code{\link{ctFit}}, from any
#'   backend.
#' @return A long-format data.frame, as originally supplied to
#'   \code{\link{ctFit}}: one row per observation, with the subject id, time
#'   and manifest/predictor columns named as in the original data.
#' @export
#' @examples
#' d <- ctFitData(ctstantestfit)
#' head(d)
ctFitData <- function(fit) {
  if (!inherits(fit, c('ctStanFit', 'ctJuliaFit'))) stop('fit object is not a ctsem fit!')
  .ctFitLongData(fit)
}
