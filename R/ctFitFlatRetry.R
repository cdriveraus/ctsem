#' Fit again from the other end of a flat parameter's range
#'
#' A fit can leave a parameter far out where the likelihood hardly changes with
#' it, such as an auto-effect so fast or so slow that its value no longer
#' matters, or a variance or random effect standard deviation at zero. A
#' higher maximum may then have that parameter elsewhere and the rest of the
#' model arranged around it. \code{\link{ctFit}} names such parameters in a
#' message. This function fits the model again from the default starting
#' values with each of them held at the other end of its range, then lets it
#' go, and keeps the better fit.
#'
#' @param fit A julia fit from \code{\link{ctFit}}.
#' @param cores Number of CPU cores to use. NULL uses the fit's own setting.
#' @param ... Other arguments passed to \code{\link{ctFit}}, as for
#'   \code{\link{ctFitUpdate}}.
#' @return The new fit where fitting again raised the log likelihood, otherwise
#'   \code{fit} unchanged.
#' @details Each parameter costs up to two more optimisations.
#'   \code{optimcontrol = list(flatretry = TRUE)} in \code{\link{ctFit}} does
#'   the same during the fit.
#' @examples
#' \dontrun{
#' fit <- ctFit(datalong, model)
#' fit <- ctFitFlatRetry(fit)
#' }
#' @export
ctFitFlatRetry <- function(fit, cores = NULL, ...) {
  if (!.ctFitIsJulia(fit)) stop("ctFitFlatRetry() applies to julia fits.",
    call. = FALSE)
  names <- as.character(fit$optim$flatretry_parameters)
  if (!length(names)) {
    message("No parameter of this fit ran off to where the data say little ",
      "about it, so the fit is returned unchanged.")
    return(fit)
  }
  dots <- list(...)
  # The fit's own call with the escape on, from its estimate (`ctFitUpdate()`):
  # the stage there ends at once, and its certification loop runs the escape
  # from the start the fit takes without `inits` (`.ctJuliaOptimiseFit()`).
  oc <- .ctJuliaOr(.ctFitReplayArgs(fit)$optimcontrol, list())
  if (!is.null(dots$optimcontrol)) oc <- utils::modifyList(oc, dots$optimcontrol)
  oc$flatretry <- TRUE
  dots$optimcontrol <- oc
  if (!is.null(cores)) dots$cores <- cores
  # The refit says again what the fit said; where nothing is kept, the fit is
  # returned as it was and so is what it said.
  said <- list()
  newfit <- withCallingHandlers(do.call(ctFitUpdate, c(list(fit, refit = TRUE),
    dots)),
    message = function(m) {
      said[[length(said) + 1L]] <<- m
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      said[[length(said) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  kept <- any(vapply(.ctJuliaOr(newfit$optim$corrections, list()), function(r)
    identical(r$status, "runoff") && isTRUE(r$resumed), logical(1)))
  if (kept) {
    for (s in said) if (inherits(s, "warning")) warning(s) else message(s)
    return(newfit)
  }
  message("Fitting again from the other end of the range of ",
    paste(names, collapse = ", "), " did not raise the log likelihood, so the ",
    "fit is returned unchanged.")
  fit
}
