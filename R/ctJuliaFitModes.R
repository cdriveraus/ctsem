# A fit's modes: one objective for everything that reads the fit.
#
# Within an optimisation each unit's inner solve starts from its mode at the
# last accepted point (the anchors, CLAUDE.md "Within an optimisation the inner
# solve is anchored"), and the run searches again where it stops and ends, so
# the estimate, its certification and its Hessian are the anchored objective's.
# What reads the fit afterwards used to clear the anchors and search every unit
# from the origin at every evaluation: a different function wherever a unit has
# more than one maximum, and on the affect pilot twice the cost (a gradient
# 115 s searched against 56 anchored, an interval check of 17 hours).
#
# So the fit keeps the modes its run ended with (`fit$laplace$modes`: units, a
# NaN-padded matrix of modes, and `at`, the free parameters they were found
# at), and every stage that evaluates the fit's objective holds them for its
# duration (`.ctJuliaHoldModes()`): each unit's Newton starts from them and
# follows each maximum from `at` to the point, the anchors never moved, so the
# value is a function of the parameters given the fit -- a density a sampler
# can target, which a chain's own history would not give wherever a unit has
# two maxima. A search from the origin is then the check, not the reference,
# and it can only add maxima: the stages whose draws reach away from the
# estimate compare the two at a few of them (`.ctJuliaModesCheck()`) and say
# where the search finds what the fit's modes miss.
#
# Held: the uncertainty stage (the Hessian where it is formed again, the
# flat-direction screen, the interval check, importance sampling and DAIS),
# the random-effects check and the refined fit's placement, the effect
# information, the certification of a stage that ended without a finish, and
# the marginal samplers. A run inside a stage sets its own anchors and clears
# them at its end. A fit without modes -- not on the Laplace route, from
# before they were kept, or every unit's inner problem concave -- is read by
# the search, as before; `ctsem_set_mode_anchors!(false)` is the oracle.

# The modes an engine run ended with (`fit_mode_*`, `_ctsem_fit_modes`), at
# its minimizer; NULL where it anchored none. `[[`: an engine result is a list
# whose names `$` would match by prefix.
.ctJuliaResultModes <- function(result) {
  units <- suppressWarnings(as.integer(result[["fit_mode_units"]]))
  values <- result[["fit_mode_values"]]
  at <- suppressWarnings(as.numeric(result[["minimizer"]]))
  if (!length(units) || !any(units >= 1L, na.rm = TRUE) || is.null(values) ||
      !length(at) || !all(is.finite(at))) return(NULL)
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  if (ncol(values) != length(units)) return(NULL)
  list(units = units, values = values, at = at)
}

# The fit's modes, or NULL.
.ctJuliaFitModes <- function(fit) {
  if (!is.list(fit$laplace)) return(NULL)
  fit$laplace[["modes"]]
}

# Holds `modes` (`.ctJuliaResultModes()`) on the objective built from `spec`
# (the hybrid, or the layer with `hybrid = FALSE`, whose units take only modes
# of their own length) until the function that called this returns. Returns
# the units held, 0 where nothing was: no modes, not the Laplace route, or the
# engine refused them, in which case the stage searches as before.
.ctJuliaHoldModes <- function(spec, modes, hybrid = TRUE, envir = parent.frame()) {
  if (is.null(modes) || is.null(spec$laplace) || isTRUE(spec$laplace$quadratic)) {
    return(invisible(0L))
  }
  handle <- if (inherits(spec, "ctJuliaModel")) spec else
    structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- .ctJuliaObjective(handle, hybrid = hybrid)
  held <- tryCatch(as.integer(.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_laplace_hold_modes!", objective,
    .ctJuliaNumericVector(modes$at), .ctJuliaVector(as.integer(modes$units)),
    modes$values))[1L], error = function(e) 0L)
  if (isTRUE(held > 0L)) {
    release <- function() try(.ctJuliaCall(
      "ContinuousTimeSEM.ctsem_laplace_release_modes!", objective, .defer = TRUE),
      silent = TRUE)
    do.call(on.exit, list(as.call(list(release)), add = TRUE), envir = envir)
  }
  invisible(if (isTRUE(held > 0L)) held else 0L)
}

# The fit's modes against a search from the origin at the rows of `points`
# (free parameters): the engine's `ctsem_laplace_modes_check`, under a hold of
# its own. NULL where no modes could be held. The record says how many points
# were checked and evaluated, at how many the search found maxima the modes
# miss worth more than `tol` nats to some unit, the most the held density
# missed at a point (`gain`, the sum over units), the units named and the most
# each missed, and the modes the search found for them (`found`, as known
# modes are kept).
.ctJuliaModesCheck <- function(spec, modes, points, hybrid = TRUE, tol = 0.1) {
  points <- as.matrix(points)
  if (!nrow(points) || .ctJuliaHoldModes(spec, modes, hybrid = hybrid) == 0L) {
    return(NULL)
  }
  handle <- if (inherits(spec, "ctJuliaModel")) spec else
    structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  columns <- t(points)
  storage.mode(columns) <- "double"
  r <- tryCatch(.ctJuliaGet(module$ctsem_laplace_modes_check(
    .ctJuliaObjective(handle, hybrid = hybrid), columns, tol = as.numeric(tol))),
    error = function(e) NULL)
  if (is.null(r)) return(NULL)
  gain <- as.numeric(r$gain)
  units <- as.integer(r$units)
  named <- units >= 1L
  found <- as.integer(r$found_units)
  list(checked = nrow(points), evaluated = sum(is.finite(gain)),
    missed = sum(is.finite(gain) & gain > tol), gain = gain,
    max_gain = if (any(is.finite(gain))) max(gain, na.rm = TRUE) else NA_real_,
    units = units[named], unit_gain = as.numeric(r$unit_gain)[named],
    found = if (any(found >= 1L)) list(units = found,
      values = as.matrix(r$found_values)) else NULL,
    tol = tol)
}

# Rows of draws from the normal at `est` with covariance `cov`, from a stream
# of their own: where an uncertainty's intervals and its interval check read
# the fit's objective.
.ctJuliaModesCheckPoints <- function(est, cov, n = 16L, seed = 20261104L) {
  A <- tryCatch(.ctCovRoot(cov), error = function(e) NULL)
  if (is.null(A) || !all(is.finite(A))) return(NULL)
  restore <- .ctStreamKeep(seed)
  on.exit(restore(), add = TRUE)
  Z <- matrix(stats::rnorm(n * ncol(A)), n, ncol(A))
  sweep(Z %*% t(A), 2, as.numeric(est), "+")
}

# One sentence where a check found maxima the fit's modes miss: what read the
# fit's objective (`what`), where, and how much.
.ctJuliaModesCheckWarn <- function(check, what) {
  if (is.null(check) || !isTRUE(check$missed > 0L)) return(invisible(NULL))
  warning(sprintf(paste0("At %d of %d %s checked, a search finds maxima of %d ",
    "unit(s) that the fit's modes miss, worth up to %.2g nats to the density ",
    "(units %s); %s read the objective at the fit's modes, which leaves them ",
    "out. Sampling the joint, ctFitUncertainty(fit, 'sample', intoverpop = ",
    "'none'), integrates over every maximum."), check$missed, check$evaluated,
    if (identical(what, "the sample")) "posterior draws" else "draws",
    length(check$units), check$max_gain,
    paste(utils::head(check$units, 5), collapse = ", "), what), call. = FALSE)
  invisible(NULL)
}
