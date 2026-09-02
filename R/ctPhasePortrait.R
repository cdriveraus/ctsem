# Phase portrait ---------------------------------------------------------------
#
# What a nonlinear continuous time model does is a property of its *vector
# field*: at each point of the state space, which way and how fast the system
# moves. Every other summary in the package reports that field through a single
# number per matrix cell, at one point, which for a linear model loses nothing
# and for a nonlinear one loses the thing that made it worth fitting.
#
# The portrait draws the field directly over a plane of the state space, with
# the structure that organises it: the nullclines, where one process stops
# changing, and their intersections, where the whole system does.
#
# Linear models are supported and drawn correctly -- the field is affine, the
# nullclines are straight, and there is exactly one fixed point -- but a linear
# field can be read off the DRIFT matrix, so the picture is not news for one.
# That is also why a stan fit can be drawn at all: a constant field needs no
# engine.

# The deterministic derivative as a function of state, however the backend can
# supply it.
#
# A julia fit re-materialises the matrices at each point, which is the only way
# to get a state-dependent field. A stan fit cannot, but a *linear* stan fit
# does not need to: its field is DRIFT x + CINT with both constant, so it comes
# from the summary matrices. A nonlinear stan fit is refused rather than drawn
# as though it were linear, which would be the one genuinely misleading outcome.
.ctFieldFunction <- function(fit, tipreds = NULL) {
  nlatent <- .ctFitNlatent(fit)
  if (inherits(fit, "ctJuliaFit")) {
    # Discrete time as below: the engine returns DRIFT x + CINT, which is the
    # next state rather than the change, so subtract the current one.
    discrete <- !isTRUE(.ctFitModelObject(fit)$continuoustime)
    return(function(state) {
      value <- .ctNonlinearDerivative(fit, state, tipreds, nlatent)$derivative
      if (discrete) value <- value - state[seq_len(nlatent)]
      value
    })
  }
  if (isTRUE(ctModelIsNonlinear(fit))) stop(call. = FALSE, paste0(
    "The vector field of a model whose matrices depend on the state can only ",
    "be evaluated by the engine that fitted it; this needs backend='julia'."))
  linear <- .ctPhaseLinearParts(fit)
  function(state) as.numeric(linear$drift %*% state[seq_len(nlatent)]) + linear$cint
}

# DRIFT and CINT as the field uses them.
#
# In discrete time ctsem's DRIFT is the one step transition, so the change per
# step is (DRIFT - I) x + CINT, not DRIFT x + CINT -- the latter is where the
# system goes next, not which way it is moving, and plotting it as an arrow
# would point everything at the origin's image rather than along the flow.
.ctPhaseLinearParts <- function(fit) {
  nlatent <- .ctFitNlatent(fit)
  mats <- suppressMessages(ctSummaryMatrices(fit))
  drift <- as.matrix(mats$DRIFT)[seq_len(nlatent), seq_len(nlatent), drop = FALSE]
  if (!isTRUE(.ctFitModelObject(fit)$continuoustime)) drift <- drift - diag(nlatent)
  list(drift = drift, cint = as.numeric(mats$CINT[seq_len(nlatent), 1]))
}

# Where to draw. 'data' follows the smoothed latent states, which is where the
# model was actually informed; 'sd' follows the model's own stationary spread,
# which is defined even when there is little data to speak of.
.ctPhaseRange <- function(fit, latents, extent, padding = 0.1) {
  if (is.matrix(extent)) {
    return(lapply(seq_len(ncol(extent)), function(i) extent[, i]))
  }
  if (is.list(extent)) return(extent)
  if (identical(extent, "sd")) {
    mats <- suppressMessages(ctSummaryMatrices(fit))
    centre <- as.numeric(mats$T0MEANS)[latents]
    spread <- sqrt(pmax(diag(as.matrix(mats$asymDIFFUSIONcov))[latents], 0))
    return(lapply(seq_along(latents), function(i)
      c(centre[i] - 2 * spread[i], centre[i] + 2 * spread[i])))
  }
  smoothed <- suppressMessages(ctKalmanArray(fit, pointest = TRUE))$etasmooth
  lapply(latents, function(k) {
    values <- as.numeric(smoothed[, , k])
    values <- values[is.finite(values)]
    span <- diff(range(values))
    range(values) + c(-1, 1) * padding * max(span, 1e-8)
  })
}

#' Phase portrait of a fitted continuous time model
#'
#' @description
#' Draw the model's vector field over a plane of the latent state space: at
#' each point, the direction and speed the deterministic part of the system
#' moves. For a model whose matrices depend on the latent state -- the reason
#' the julia backend exists -- this shows behaviour that no single reported
#' DRIFT matrix can, because there is no single DRIFT matrix.
#'
#' @details
#' Three things are drawn over the field, each of which can be turned off:
#'
#' \itemize{
#'   \item \strong{Nullclines}, the curves where one process stops changing.
#'     Where the nullclines of two processes cross, the whole system stops:
#'     that is a fixed point, and a nonlinear system may have several.
#'   \item \strong{The fixed point} found by \code{state='asymptotic'}, marked
#'     when it lies inside the plotted region. Only one is found, by Newton from
#'     T0MEANS; the nullcline crossings show whether there are others.
#'   \item \strong{Smoothed trajectories} for a few subjects, so the field can
#'     be read against where the data actually went.
#' }
#'
#' With more than two latent processes the remaining ones are held at
#' \code{state}, so the picture is a slice rather than the whole story --
#' change \code{state} and the slice changes.
#'
#' A linear model has an affine field, straight nullclines and exactly one
#' fixed point. It is drawn correctly, from a stan fit as readily as a julia
#' one, but it shows nothing the DRIFT matrix did not already say. A nonlinear
#' stan fit is refused rather than drawn as if it were linear.
#'
#' @param fit A ctJuliaFit, or a ctStanFit of a linear model.
#' @param latents Two latent processes, by name or index, forming the plane.
#' @param extent \code{'data'} (the default) spans the smoothed latent states,
#'   \code{'sd'} spans two stationary standard deviations either side of
#'   T0MEANS, or supply a list of two numeric ranges.
#' @param gridsize Arrows per side. The field costs one engine call per point,
#'   so cost is quadratic in this.
#' @param state Where the *other* latent processes are held, as in
#'   \code{\link{ctDiscretePars}}: \code{'asymptotic'}, \code{'mean'},
#'   \code{'T0MEANS'}, or a numeric state.
#' @param tipreds Time independent predictor values, or NULL for the population.
#' @param nullclines Draw the curves where each process stops changing?
#' @param fixedpoint Mark the asymptotic state when it is inside the region?
#' @param trajectories Number of subjects' smoothed trajectories to overlay, or
#'   0 for none.
#' @param plot If FALSE, return the field, nullclines and trajectories as data
#'   rather than a plot.
#' @param ... Ignored.
#' @return A ggplot, or a list of data frames when \code{plot=FALSE}.
#' @examples
#' # A linear stan fit needs no engine, so the bundled example fit is enough.
#' ctPhasePortrait(ctstantestfit, gridsize = 7, plot = FALSE)
#' @seealso \code{\link{ctContextDependence}}, \code{\link{ctDiscretePars}}
#' @export
ctPhasePortrait <- function(fit, latents = 1:2, extent = "data", gridsize = 15,
  state = "asymptotic", tipreds = NULL, nullclines = TRUE, fixedpoint = TRUE,
  trajectories = 0, plot = TRUE, ...) {

  model <- .ctFitModelObject(fit)
  latentNames <- model$latentNames
  if (is.character(latents)) latents <- match(latents, latentNames)
  latents <- as.integer(latents)
  if (length(latents) != 2L || any(is.na(latents))) {
    stop("latents must name or index exactly two latent processes.", call. = FALSE)
  }
  nlatent <- .ctFitNlatent(fit)

  # The plane is drawn through a point; everything not on the plane stays there.
  base <- if (inherits(fit, "ctJuliaFit")) {
    resolved <- try(.ctResolveState(fit, state, tipreds = tipreds), silent = TRUE)
    if (inherits(resolved, "try-error") || is.null(resolved$state)) {
      .ctContextBaseState(fit, tipreds)
    } else resolved$state
  } else {
    as.numeric(suppressMessages(ctSummaryMatrices(fit))$T0MEANS)[seq_len(nlatent)]
  }

  field <- .ctFieldFunction(fit, tipreds)
  ranges <- .ctPhaseRange(fit, latents, extent)
  xs <- seq(ranges[[1L]][1L], ranges[[1L]][2L], length.out = gridsize)
  ys <- seq(ranges[[2L]][1L], ranges[[2L]][2L], length.out = gridsize)

  grid <- expand.grid(x = xs, y = ys)
  derivatives <- matrix(NA_real_, nrow(grid), 2L)
  for (i in seq_len(nrow(grid))) {
    point <- base
    point[latents[1L]] <- grid$x[i]
    point[latents[2L]] <- grid$y[i]
    value <- try(field(point), silent = TRUE)
    if (!inherits(value, "try-error")) derivatives[i, ] <- value[latents]
  }
  grid$dx <- derivatives[, 1L]
  grid$dy <- derivatives[, 2L]
  grid$speed <- sqrt(grid$dx^2 + grid$dy^2)

  # Arrows are scaled to the grid, not to the data: what matters visually is
  # direction and relative speed, and an unscaled field either vanishes or
  # covers the plot depending on the units the processes happen to be in.
  cell <- min(diff(xs)[1L], diff(ys)[1L])
  longest <- max(grid$speed, na.rm = TRUE)
  scaling <- if (is.finite(longest) && longest > 0) 0.9 * cell / longest else 0
  grid$xend <- grid$x + grid$dx * scaling
  grid$yend <- grid$y + grid$dy * scaling

  out <- list(field = grid, latents = latentNames[latents])

  if (isTRUE(nullclines)) {
    # A nullcline is a zero contour of one component of the field, so the
    # contouring routine finds it without the field being solved anywhere.
    out$nullclines <- do.call(rbind, lapply(seq_len(2L), function(component) {
      z <- matrix(if (component == 1L) grid$dx else grid$dy, gridsize, gridsize)
      if (any(!is.finite(z))) return(NULL)
      lines <- try(grDevices::contourLines(xs, ys, z, levels = 0), silent = TRUE)
      if (inherits(lines, "try-error") || !length(lines)) return(NULL)
      do.call(rbind, lapply(seq_along(lines), function(k) data.frame(
        x = lines[[k]]$x, y = lines[[k]]$y,
        process = latentNames[latents[component]],
        piece = paste0(component, "_", k), stringsAsFactors = FALSE)))
    }))
  }

  if (isTRUE(fixedpoint)) {
    point <- if (inherits(fit, "ctJuliaFit")) {
      try(.ctContextAsymptoticState(fit, tipreds = tipreds), silent = TRUE)
    } else {
      # A linear field has its fixed point in closed form. Solved from the same
      # DRIFT and CINT the field uses, not read from asymCINT: the summary
      # collapses each matrix over draws separately, so the reported asymCINT
      # is the median of the asymptotes rather than the asymptote of the median
      # matrices, and does not sit where the drawn field vanishes.
      try({
        linear <- .ctPhaseLinearParts(fit)
        as.numeric(solve(-linear$drift, linear$cint))
      }, silent = TRUE)
    }
    if (!inherits(point, "try-error") && length(point) >= max(latents) &&
      all(is.finite(point[latents]))) {
      inside <- point[latents[1L]] >= min(xs) && point[latents[1L]] <= max(xs) &&
        point[latents[2L]] >= min(ys) && point[latents[2L]] <= max(ys)
      if (inside) out$fixedpoint <- data.frame(x = point[latents[1L]],
        y = point[latents[2L]])
    }
  }

  if (trajectories > 0) {
    kalman <- suppressMessages(ctKalmanArray(fit, pointest = TRUE))
    eta <- kalman$etasmooth
    subject <- if (!is.null(kalman$id)) as.integer(kalman$id) else
      rep(1L, dim(eta)[2L])
    keep <- unique(subject)[seq_len(min(trajectories, length(unique(subject))))]
    rows <- which(subject %in% keep)
    times <- if (!is.null(kalman$time)) as.numeric(kalman$time)[rows] else seq_along(rows)
    # Ordered by time within subject: geom_path joins points in row order, and
    # a path through a phase plane in the wrong order is not a trajectory.
    rows <- rows[order(subject[rows], times)]
    out$trajectories <- data.frame(
      x = as.numeric(eta[1L, rows, latents[1L]]),
      y = as.numeric(eta[1L, rows, latents[2L]]),
      Subject = factor(subject[rows]))
  }

  if (!isTRUE(plot)) return(out)
  .ctPhasePortraitPlot(out, fit)
}

.ctPhasePortraitPlot <- function(portrait, fit) {
  x <- y <- xend <- yend <- speed <- process <- piece <- Subject <- NULL
  contextual <- isTRUE(ctModelIsNonlinear(fit))

  g <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = portrait$field,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = speed),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.14, "cm")), linewidth = 0.35) +
    ggplot2::scale_colour_viridis_c(name = "Speed")

  if (!is.null(portrait$nullclines) && nrow(portrait$nullclines)) {
    g <- g + ggplot2::geom_path(data = portrait$nullclines,
      ggplot2::aes(x = x, y = y, group = piece, linetype = process),
      linewidth = 0.6, colour = "grey25") +
      ggplot2::labs(linetype = "Nullcline of")
  }
  if (!is.null(portrait$trajectories) && nrow(portrait$trajectories)) {
    g <- g + ggplot2::geom_path(data = portrait$trajectories,
      ggplot2::aes(x = x, y = y, group = Subject), colour = "grey45",
      alpha = 0.5, linewidth = 0.35)
  }
  if (!is.null(portrait$fixedpoint)) {
    g <- g + ggplot2::geom_point(data = portrait$fixedpoint,
      ggplot2::aes(x = x, y = y), shape = 21, size = 3, stroke = 1,
      fill = "white", colour = "black")
  }

  g + ggplot2::labs(
    x = portrait$latents[1L], y = portrait$latents[2L],
    title = if (contextual) "Phase portrait (state dependent dynamics)" else
      "Phase portrait (linear dynamics)",
    subtitle = paste0("Arrows: deterministic change. Lines: nullclines",
      if (!is.null(portrait$fixedpoint)) ". Circle: fixed point" else "", ".")) +
    ggplot2::theme_minimal()
}
