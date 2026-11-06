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
# Held: the uncertainty stage (the Hessian where it is formed again, the
# flat-direction screen, the interval check, importance sampling and DAIS),
# the random-effects check and the refined fit's placement, the effect
# information, the certification of a stage that ended without a finish, and
# the marginal samplers. One function for all of them, the one the estimate
# maximises: a Hessian formed by differences cannot see a unit land on
# different maxima either side of the estimate, and a check centred on the
# estimate reads the shape of that function's posterior rather than the gap
# between two functions. It is also the cheapest: on the affect pilot a
# gradient at draws a standard error out took 81-85 s held, 117-129 s searched
# from the origin (16 threads).
#
# Away from the estimate a search finds maxima the held modes miss, and misses
# some they keep: at three draws from the pilot's normal the two differed by
# 10.0, -8.5 and 4.7 nats, each finding up to 25 nats the other did not.
# Neither is the integral. Combining them was tried and dropped (2026-11-06):
# a maximum added beside the best, across a shallow valley, fades the best's
# weight in the mode sum (`_laplace_mode_sum`), so the combination fell below
# the held value at two of the three draws, at twice the held cost. Taking
# each unit's best of the two instead jumps wherever a search from the origin
# changes basin. A unit whose term depends on which maxima are found is one
# Laplace serves poorly, and sampling the joint, the default for intervals to
# report, finds no modes at all.
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
# of their own length) until the function that called this returns. A hold
# inside another lasts until its own caller returns, and the one under it then
# holds again. Returns the units held, 0 where nothing was: no modes, not the
# Laplace route, or the engine refused them, in which case the stage searches
# as before.
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
