# A fit's modes: one objective for everything that reads the fit.
#
# Within an optimisation each unit's inner solve starts from its mode at the
# last accepted point (the anchors, CLAUDE.md "Within an optimisation the inner
# solve is anchored"), and the run searches again where it stops and ends, so
# the estimate, its certification and its Hessian are the anchored objective's.
# What reads the fit afterwards used to clear the anchors and search every unit
# from the origin at every evaluation: a different function wherever a unit has
# more than one maximum.
#
# So the fit keeps the modes its run ended with (`fit$laplace$modes`: units, a
# NaN-padded matrix of modes, and `at`, the free parameters they were found
# at), and every stage that evaluates the fit's objective holds them for its
# duration (`.ctJuliaHoldModes()`): each unit's Newton starts from them and
# follows each maximum from `at` to the point, the anchors never moved, so the
# value is a function of the parameters given the fit -- a density a sampler
# can target, which a chain's own history would not give wherever a unit has
# two maxima.
#
# Near the estimate the followed modes are the run's function, a gradient at
# about two thirds of a search's cost on the affect pilot. Further out neither
# finds every maximum: at three draws from the pilot's normal the followed
# modes and a search from the origin differed by 10.0, -8.5 and 4.7 nats, and a
# search found maxima the modes missed worth up to 25 nats at each draw. So the
# stages whose points reach far hold with `search = TRUE`, and each unit's term
# is the combination over the maxima of both (`_laplace_mode_sum`), the term as
# the objective defines it over every maximum either finds. That is not always
# above both: a maximum the other set adds beside the best, across a shallow
# valley, fades the best's weight, and at those three draws the combination
# came out 2.55 below the followed value, 11.3 above both, and 1.13 below the
# followed value. The search-and-follow gradient there took 161-175 s against
# 117-129 s searched and 81-85 s followed (16 threads). The stages: the
# uncertainty stage after
# its Hessian (the flat-direction screen, the interval check, importance
# sampling and DAIS), the effect information, whose passes set the other
# effects' loadings to zero, and the marginal samplers. Those that stay at the
# estimate follow alone: the Hessian where formed again, the random-effects
# check and placement, and the certification of a stage that ended without a
# finish.
#
# A run inside a stage sets its own anchors and clears them at its end. A fit
# without modes -- not on the Laplace route, from before they were kept, or
# every unit's inner problem concave -- is read by the search, as before;
# `ctsem_set_mode_anchors!(false)` is the oracle.

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
# of their own length) until the function that called this returns; with
# `search`, each unit also searches from the origin and combines the maxima of
# both. A hold inside another lasts until its own caller returns, and the one
# under it then holds again as it was. Returns the units held, 0 where nothing
# was: no modes, not the Laplace route, or the engine refused them, in which
# case the stage searches as before.
.ctJuliaHoldModes <- function(spec, modes, hybrid = TRUE, search = FALSE,
  envir = parent.frame()) {
  if (is.null(modes) || is.null(spec$laplace) || isTRUE(spec$laplace$quadratic)) {
    return(invisible(0L))
  }
  handle <- if (inherits(spec, "ctJuliaModel")) spec else
    structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- .ctJuliaObjective(handle, hybrid = hybrid)
  held <- tryCatch(as.integer(.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_laplace_hold_modes!", objective,
    .ctJuliaNumericVector(modes$at), .ctJuliaVector(as.integer(modes$units)),
    modes$values, search = isTRUE(search)))[1L], error = function(e) 0L)
  if (isTRUE(held > 0L)) {
    release <- function() try(.ctJuliaCall(
      "ContinuousTimeSEM.ctsem_laplace_release_modes!", objective, .defer = TRUE),
      silent = TRUE)
    do.call(on.exit, list(as.call(list(release)), add = TRUE), envir = envir)
  }
  invisible(if (isTRUE(held > 0L)) held else 0L)
}
