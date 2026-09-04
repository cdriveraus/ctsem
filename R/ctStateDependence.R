# How a parameter varies with the state ----------------------------------------
#
# The most direct picture of a state-dependent model, and the one every other
# summary is a single sample from. `summary()` reports DRIFT[1,1] as one number
# at one state; this plots it across the range the process actually occupies,
# with the posterior spread at each point.
#
# It is cheap, because the evaluation point is the only thing that varies along
# the grid: the whole posterior can be materialised at one state in a single
# engine call, so the cost is one call per grid point rather than one per point
# per draw.

# Which matrices are worth plotting: the model matrices a reader recognises,
# excluding the Jacobian blocks, which are derivatives rather than parameters.
.ctStateDependenceMatrices <- c("DRIFT", "CINT", "LAMBDA", "MANIFESTMEANS",
  "DIFFUSIONcov", "MANIFESTcov", "TDPREDEFFECT", "T0MEANS", "PARS")

#' How model matrix cells vary with a latent state
#'
#' @description
#' Plot every context-dependent matrix cell against one latent process, across
#' the range that process occupies, with the posterior spread at each point.
#'
#' This is the plot that says what a state-dependent model actually claims.
#' \code{summary()} reports each cell as a single number at a single state;
#' here the same cell is drawn as the function of the state it really is, and
#' the summary's number is one point on that curve.
#'
#' @details
#' Cells are chosen by measuring rather than by parsing: every model matrix is
#' materialised across the grid, and a cell is kept when its value actually
#' moves. That catches dependence arriving indirectly -- through a PARS entry,
#' or through a covariance built from a state-dependent cell -- without having
#' to trace it symbolically.
#'
#' A model with no state-dependent cells has nothing to draw, and says so
#' rather than returning an empty plot.
#'
#' @param fit A ctJuliaFit.
#' @param along The latent process to vary, by name or index.
#' @param extent \code{'data'} (the default) spans the smoothed latent states,
#'   \code{'sd'} spans two stationary standard deviations either side of
#'   T0MEANS, or supply a numeric range of length two.
#' @param gridsize Points along the range. One engine call each.
#' @param nsamples Posterior draws for the interval, or 1 for the point
#'   estimate alone.
#' @param probs Three quantiles: lower ribbon, line, upper ribbon.
#' @param state Where the other latent processes are held: \code{'asymptotic'},
#'   \code{'mean'}, \code{'T0MEANS'}, or a numeric state.
#' @param tipreds Time independent predictor values, or NULL for the population.
#' @param tolerance A cell is treated as varying when its range across the grid
#'   exceeds this fraction of its own magnitude.
#' @param plot If FALSE, return the values as a data frame rather than a plot.
#' @return A ggplot, or a data frame when \code{plot=FALSE}.
#' @examples
#' \donttest{
#' # Needs a working julia backend, so the example is inert where Julia is
#' # absent -- including on CRAN, whose check machines have none.
#' if (isTRUE(ctJuliaStatus()$available)) {
#'   set.seed(1)
#'   generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
#'     manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
#'     LAMBDA = diag(2), DRIFT = matrix(c(-.4, .1, 0, -.3), 2, 2),
#'     CINT = matrix(c(.2, .1), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
#'     MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c(.5, 0, 0, .4), 2, 2)))
#'   datalong <- as.data.frame(suppressMessages(ctGenerate(generating, n.subjects = 25,
#'     burnin = 5, dtmean = 1, logdtsd = .1, wide = FALSE, Tpoints = 12)))
#'
#'   # eta1's own decay depends on where eta2 is -- a state-dependent DRIFT cell.
#'   model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
#'     manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
#'     LAMBDA = diag(2), PARS = c('dr11|-log1p_exp(param)'),
#'     DRIFT = matrix(c('dr11 * (1 + 0.2 * eta2)', 'd21', 0, 'd22'), 2, 2),
#'     CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
#'     MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
#'   model$pars$indvarying <- FALSE
#'
#'   fit <- ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
#'   ctStateDependencePlot(fit, along = 'eta2', gridsize = 11, nsamples = 5)
#' }
#' }
#' @seealso \code{\link{ctContextDependence}}, \code{\link{ctPhasePortrait}}
#' @export
ctStateDependencePlot <- function(fit, along = 1, extent = "data", gridsize = 25,
  nsamples = 20, probs = c(.025, .5, .975), state = "asymptotic", tipreds = NULL,
  tolerance = 1e-6, plot = TRUE) {

  if (!inherits(fit, "ctJuliaFit")) stop(call. = FALSE, paste0(
    "Evaluating a cell across a range of states needs the engine that fitted ",
    "the model; this requires backend='julia'."))

  model <- .ctFitModelObject(fit)
  latentNames <- model$latentNames
  if (is.character(along)) along <- match(along, latentNames)
  along <- as.integer(along)
  if (length(along) != 1L || is.na(along)) {
    stop("along must name or index one latent process.", call. = FALSE)
  }
  if (!nrow(.ctFitConditionalCells(fit))) {
    stop(call. = FALSE, paste0("No cell of this model depends on the latent ",
      "state or a time dependent predictor, so there is nothing to plot. ",
      "ctModelIsNonlinear() says the same, without the error."))
  }

  resolved <- .ctResolveState(fit, state, tipreds = tipreds)
  base <- if (is.null(resolved$state)) .ctContextBaseState(fit, tipreds) else resolved$state
  span <- if (is.numeric(extent) && length(extent) == 2L) extent else
    .ctPhaseRange(fit, along, extent)[[1L]]
  grid <- seq(span[1L], span[2L], length.out = gridsize)

  draws <- .ctBackendRawSamples(fit)
  nsamples <- max(1L, min(as.integer(nsamples), nrow(draws)))
  if (nsamples < nrow(draws)) {
    draws <- draws[round(seq(1, nrow(draws), length.out = nsamples)), , drop = FALSE]
  }

  # One engine call per grid point covers the whole posterior at that point:
  # `state` is a keyword of the materialisation, the draws are its columns.
  collected <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    point <- base
    point[along] <- grid[i]
    arrays <- try(.ctBackendPopArrays(fit, samples = draws, state = point,
      tipreds = tipreds), silent = TRUE)
    if (inherits(arrays, "try-error")) next
    collected[[i]] <- arrays
  }
  keep <- !vapply(collected, is.null, logical(1L))
  if (!any(keep)) stop("The engine could not materialise the matrices anywhere on this range.",
    call. = FALSE)
  grid <- grid[keep]
  collected <- collected[keep]

  values <- do.call(rbind, lapply(seq_along(grid), function(i) {
    arrays <- collected[[i]]
    do.call(rbind, lapply(.ctStateDependenceMatrices, function(name) {
      block <- arrays[[paste0("pop_", name)]]
      if (is.null(block) || length(dim(block)) != 3L) return(NULL)
      quantiles <- apply(block, c(2L, 3L), stats::quantile, probs = probs,
        na.rm = TRUE)
      expand <- expand.grid(row = seq_len(dim(block)[2L]), col = seq_len(dim(block)[3L]))
      data.frame(matrix = name, row = expand$row, col = expand$col,
        along = grid[i],
        lower = as.numeric(quantiles[1L, , ]), middle = as.numeric(quantiles[2L, , ]),
        upper = as.numeric(quantiles[3L, , ]), stringsAsFactors = FALSE)
    }))
  }))

  # Keep the cells that move. Measured rather than parsed, so dependence
  # arriving through a PARS entry or through a covariance built from a
  # state-dependent cell is caught without tracing it symbolically.
  values$cell <- paste0(values$matrix, "[", values$row, ",", values$col, "]")
  varying <- vapply(split(values$middle, values$cell), function(x) {
    scale <- max(abs(x), na.rm = TRUE)
    if (!is.finite(scale) || scale == 0) return(FALSE)
    diff(range(x, na.rm = TRUE)) / scale > tolerance
  }, logical(1L))
  values <- values[values$cell %in% names(varying)[varying], , drop = FALSE]
  if (!nrow(values)) stop(call. = FALSE, paste0(
    "No reported matrix cell changed measurably across this range of ",
    latentNames[along], ". Try along= another process, or a wider extent="))

  attr(values, "along") <- latentNames[along]
  attr(values, "stateLabel") <- resolved$label
  if (!isTRUE(plot)) return(values)

  cell <- lower <- upper <- middle <- NULL
  ggplot2::ggplot(values, ggplot2::aes(x = along)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = middle), linewidth = 0.7) +
    ggplot2::facet_wrap(~cell, scales = "free_y") +
    ggplot2::labs(x = latentNames[along], y = "Parameter value",
      title = paste0("Matrix cells as functions of ", latentNames[along]),
      subtitle = paste0("Other processes held at ", resolved$label,
        ". Ribbon: ", probs[1L] * 100, "-", probs[3L] * 100, "% of draws.")) +
    ggplot2::theme_minimal()
}
