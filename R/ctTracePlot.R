# The trace of how a fit reached its estimate, drawn.
#
# For an optimized fit, the optimiser's: what a reader wants from it is not the
# objective's final value -- the summary has that -- but its *shape*: whether
# the log posterior flattened long before the iteration budget ran out, or was
# still climbing when it stopped, and whether the gradient norm came down with
# it. Those two together are what distinguish a fit that converged from one
# that merely stopped, and they are the reason the trace is recorded at every
# iteration rather than at the reporting cadence. The gradient norm is drawn
# on a log scale because that is the only scale it is legible on: it routinely
# spans ten orders of magnitude, and on a linear axis every point after the
# third sits on the floor.
#
# For a sampled fit, each chain's draws of the population parameters, on the
# natural scale summary() reports them on: the picture that says whether the
# chains settled on one distribution, which R-hat summarises in a number.

#' Plot how a fit reached its estimate
#'
#' For an optimized fit, draws the log posterior and the gradient norm against
#' iteration, from the record the engine kept while fitting. For a sampled fit,
#' draws each chain's draws of every population mean and standard deviation,
#' on the scale \code{summary()} reports them on, so chains that have not
#' settled on one distribution can be seen.
#'
#' @param fit A fit from \code{\link{ctFit}} with the Julia backend.
#' @param which Which panels to draw: trace columns (\code{'objective'},
#'   \code{'gradient_norm'}, \code{'inner_converged'}) for an optimized fit,
#'   parameter names for a sampled one. Defaults to all of them.
#' @param plot If \code{FALSE}, return the data rather than a plot: the trace
#'   for an optimized fit, one row per draw and parameter for a sampled one.
#'
#' @return A ggplot, or with \code{plot = FALSE} a data frame.
#'
#' @details The Julia backend records the optimisation trace during the fit and
#'   returns it as \code{fit$optim$trace}. For output while a fit is still
#'   running, see \code{optimcontrol$callback} in \code{\link{ctFit}}. The
#'   chains of a sampled stan fit are drawn by \code{plot(fit, types =
#'   'trace')}.
#'
#' @examples
#' \donttest{
#' # Needs a working julia backend, so the example is inert where Julia is
#' # absent -- including on CRAN, whose check machines have none.
#' if (isTRUE(ctJuliaStatus()$available)) {
#'   model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
#'     LAMBDA = matrix(1))
#'   model$pars$indvarying <- FALSE
#'   fit <- ctFit(ctstantestdat, model, backend = 'julia', cores = 1)
#'   ctTracePlot(fit)
#'   # Julia keeps running between calls, holding what it has compiled;
#'   # this ends it.
#'   ctJuliaKill("session")
#' }
#' }
#' @export
ctTracePlot <- function(fit, which = NULL, plot = TRUE) {
  if (inherits(fit, 'ctStanFit')) stop("ctTracePlot() reads the trace the ",
    "julia engine records. For a sampled stan fit, plot(fit, types = 'trace') ",
    "draws its chains.", call. = FALSE)
  if (!is.data.frame(fit) && !is.null(fit$sample))
    return(.ctTracePlotChains(fit, which, plot))
  trace <- if (is.data.frame(fit)) fit else fit$optim$trace
  if (is.null(trace) || !nrow(trace)) {
    stop("This fit carries no optimisation trace. The Julia backend records ",
      "one; the Stan backend does not.", call. = FALSE)
  }
  columns <- setdiff(names(trace), "iteration")
  if (!is.null(which)) columns <- intersect(columns, which)
  if (!length(columns)) stop("No trace columns to plot.", call. = FALSE)
  if (!isTRUE(plot)) return(trace[, c("iteration", columns), drop = FALSE])

  labels <- c(objective = "log posterior", gradient_norm = "gradient norm (log10)",
    inner_converged = "inner modes converged")
  long <- do.call(rbind, lapply(columns, function(column) {
    values <- as.numeric(trace[[column]])
    label <- if (column %in% names(labels)) labels[[column]] else column
    if (identical(column, "gradient_norm")) {
      if (all(values > 0)) values <- log10(values) else label <- "gradient norm"
    }
    data.frame(iteration = trace$iteration, value = values, panel = label)
  }))
  long$panel <- factor(long$panel, levels = unique(long$panel))
  if (1 == 99) iteration <- value <- panel <- NULL
  ggplot(long, aes(x = iteration, y = value)) +
    geom_line() +
    geom_point(data = long[!duplicated(long$panel, fromLast = TRUE), ], size = 1.5) +
    facet_wrap(vars(panel), ncol = 1, scales = 'free_y') +
    labs(x = 'Iteration', y = NULL) +
    theme_bw()
}

# Each chain's draws of the population means and sds, natural scale. The rows
# of `estimate$rawposterior` are the chains' draws one chain after another,
# which `.ctBackendSummaryChains()` recovers -- and refuses when the order no
# longer holds, as after ctParticleCorrect() reweights them.
.ctTracePlotChains <- function(fit, which, plot) {
  samples <- .ctBackendRawSamples(fit)
  chains <- .ctBackendSummaryChains(fit, samples)
  if (is.null(chains)) stop("This fit's draws no longer carry the order of the ",
    "chains they came from (ctParticleCorrect() reweights them, for instance), ",
    "so there is no chain to trace.", call. = FALSE)
  lengths <- if (length(chains) > 1L) chains else
    rep(nrow(samples) %/% chains, chains)
  draws <- list(.ctBackendPopMeanSamples(fit, samples = samples)$values)
  levels <- .ctBackendConstrained(fit, samples)$randomeffectlevels
  for (l in seq_along(levels)) {
    sds <- levels[[l]]$popsd
    if (is.null(sds) || !ncol(sds)) next
    colnames(sds) <- paste0(if (length(levels) > 1L)
      paste0('sd(', levels[[l]]$level, ') ') else 'sd ', colnames(sds))
    draws[[length(draws) + 1L]] <- sds
  }
  draws <- do.call(cbind, draws)
  if (!is.null(which)) {
    unknown <- setdiff(which, colnames(draws))
    if (length(unknown)) stop("Not a population parameter of this fit: ",
      paste(unknown, collapse = ', '), ". Choose from: ",
      paste(colnames(draws), collapse = ', '), ".", call. = FALSE)
    draws <- draws[, which, drop = FALSE]
  }
  long <- data.frame(
    draw = rep(unlist(lapply(lengths, seq_len)), ncol(draws)),
    chain = factor(rep(rep(seq_along(lengths), lengths), ncol(draws))),
    parameter = factor(rep(colnames(draws), each = nrow(draws)),
      levels = colnames(draws)),
    value = as.numeric(draws))
  if (!isTRUE(plot)) return(long)
  if (1 == 99) draw <- value <- chain <- parameter <- NULL
  ggplot(long, aes(x = draw, y = value, colour = chain)) +
    geom_line(alpha = .7, linewidth = .3) +
    facet_wrap(vars(parameter), scales = 'free_y') +
    labs(x = 'Draw (after warmup)', y = NULL, colour = 'Chain') +
    theme_bw()
}
