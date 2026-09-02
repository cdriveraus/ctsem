# The optimisation trace, drawn.
#
# What a reader wants from this is not the objective's final value -- the
# summary has that -- but its *shape*: whether the log posterior flattened long
# before the iteration budget ran out, or was still climbing when it stopped,
# and whether the gradient norm came down with it. Those two together are what
# distinguish a fit that converged from one that merely stopped, and they are
# the reason the trace is recorded at every iteration rather than at the
# reporting cadence.
#
# The gradient panel is on a log scale because that is the only scale it is
# legible on: it routinely spans ten orders of magnitude between the first
# iteration and the last, and on a linear axis every point after the third sits
# on the floor.

#' Plot a fit's optimisation trace
#'
#' Draws the log posterior and the gradient norm against iteration, from the
#' record the engine kept while fitting.
#'
#' @param fit A fit from \code{\link{ctFit}} with the Julia backend.
#' @param which Which panels to draw. Defaults to every column the trace has
#'   besides the iteration index.
#' @param ... Passed to \code{plot}.
#'
#' @return The trace, invisibly.
#'
#' @details Available on the Julia backend, which records the trace during the
#'   fit and returns it as \code{fit$trace}. For output while a fit is still
#'   running, see \code{optimcontrol$callback} in \code{\link{ctFit}}.
#'
#' @export
ctTracePlot <- function(fit, which = NULL, ...) {
  trace <- if (is.data.frame(fit)) fit else fit$trace
  if (is.null(trace) || !nrow(trace)) {
    stop("This fit carries no optimisation trace. The Julia backend records ",
      "one; the Stan backend does not.", call. = FALSE)
  }
  columns <- setdiff(names(trace), "iteration")
  if (!is.null(which)) columns <- intersect(columns, which)
  if (!length(columns)) stop("No trace columns to plot.", call. = FALSE)

  labels <- c(objective = "log posterior", gradient_norm = "gradient norm",
    inner_converged = "inner modes converged")
  old <- graphics::par(mfrow = c(length(columns), 1L),
    mar = c(4, 4.5, 1.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (column in columns) {
    values <- trace[[column]]
    # The gradient norm falls across many orders of magnitude, so on a linear
    # axis everything after the first few iterations is indistinguishable from
    # zero. Only taken when the values allow it.
    logscale <- identical(column, "gradient_norm") && all(values > 0)
    graphics::plot(trace$iteration, values, type = "l",
      log = if (logscale) "y" else "",
      xlab = "iteration",
      ylab = if (column %in% names(labels)) labels[[column]] else column, ...)
    graphics::points(trace$iteration[nrow(trace)], values[nrow(trace)],
      pch = 16, cex = 0.8)
  }
  invisible(trace)
}
