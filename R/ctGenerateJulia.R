# Generating data through the julia engine.
#
# `ctGenerate()`'s own generator integrates the linear system directly with a
# matrix exponential. That is exact for what it covers and cannot cover anything
# else: a state-dependent drift, a parameter that is a function of another, any
# of the nonlinear forms the engine filters happily. The engine already
# generates -- `ctsem_generate` draws each row from its own prior predictive as
# the filter reaches it -- so the gap is reachability, not capability.
#
# The standard normals are drawn on the R side, so `set.seed()` means what a
# user expects and the two paths agree for a given seed.

# Values for parameters that generation has to pin down.
#
# `ctModeltoNumeric()` sets every free parameter to zero, which is right for a
# location and wrong for everything else: a zero DIFFUSION is a process with no
# innovation, a zero DRIFT is a random walk rather than the mean-reverting
# process the user drew, a zero MANIFESTVAR is noiseless measurement, and a zero
# LAMBDA disconnects the manifest from the latent it loads on. Data generated
# under those is degenerate in ways that are not obvious until an analysis of it
# behaves strangely.
#
# So: zero where zero is the natural centre, and a plausible non-zero elsewhere.
# These are defaults for *simulation*, not estimates of anything, and the point
# is only that data generated from them looks like data rather than like an
# artefact. Anything a user cares about they should set.
#' @keywords internal
.ctGenerateDefaults <- function() list(
  DRIFT = list(diagonal = -0.5, offdiagonal = 0),
  DIFFUSION = list(diagonal = 1, offdiagonal = 0),
  T0VAR = list(diagonal = 1, offdiagonal = 0),
  MANIFESTVAR = list(diagonal = 0.5, offdiagonal = 0),
  # A free loading defaults to one: zero would cut the manifest off from its
  # latent entirely, which is never what someone leaving LAMBDA free meant.
  LAMBDA = list(diagonal = 1, offdiagonal = 0),
  T0MEANS = list(diagonal = 0, offdiagonal = 0),
  MANIFESTMEANS = list(diagonal = 0, offdiagonal = 0),
  CINT = list(diagonal = 0, offdiagonal = 0),
  TDPREDEFFECT = list(diagonal = 0, offdiagonal = 0),
  PARS = list(diagonal = 0, offdiagonal = 0))

#' @keywords internal
.ctGenerateResolveFree <- function(model, quiet = FALSE) {
  defaults <- .ctGenerateDefaults()
  pars <- model$pars
  free <- which(is.na(pars$value))
  if (!length(free)) return(model)
  filled <- character()
  for (i in free) {
    matrix_name <- as.character(pars$matrix[i])
    spec <- defaults[[matrix_name]]
    # An unlisted matrix gets zero, which is the old behaviour and the right
    # fallback: a matrix this function has no opinion about is one where a
    # non-zero guess would be a fabrication.
    value <- if (is.null(spec)) 0 else
      if (isTRUE(pars$row[i] == pars$col[i])) spec$diagonal else spec$offdiagonal
    pars$value[i] <- value
    if (value != 0) filled <- c(filled, sprintf("%s[%d,%d]=%s", matrix_name,
      pars$row[i], pars$col[i], format(value)))
  }
  model$pars <- pars
  if (!quiet) {
    message(length(free), " free parameter", if (length(free) > 1) "s" else "",
      " had no value and were set for generation",
      if (length(filled)) paste0(", including ",
        paste(utils::head(filled, 6), collapse = ", "),
        if (length(filled) > 6) ", ..." else "") else " to zero",
      ". Set them in the model if they matter.")
  }
  model
}

# A long-format skeleton with the requested subjects and times and no
# observations: the shape generation fills in.
#' @keywords internal
.ctGenerateSkeleton <- function(model, n.subjects, times) {
  rows <- do.call(rbind, lapply(seq_len(n.subjects), function(i)
    data.frame(id = i, time = times[[i]])))
  # Zero, not NA. The engine generates only where an observation exists --
  # `_generate_row!` is handed the observed indices and writes only those -- so
  # the skeleton's *missingness pattern* is the input and its values are not:
  # the generated draw overwrites the value and the filter's update conditions
  # on the draw, never on what was there. An all-NA skeleton therefore produces
  # an all-NA result, which is what it did before this line said zero.
  for (nm in model$manifestNames) rows[[nm]] <- 0
  # Predictors at zero. A time-dependent predictor that is never non-zero has
  # no effect on the generated data, which is the honest default for a value
  # the caller did not supply; `ctGenerate`'s own path uses TDPREDMEANS and
  # this could too once the specification carries them.
  for (nm in model$TDpredNames) rows[[nm]] <- 0
  for (nm in model$TIpredNames) rows[[nm]] <- 0
  rows
}

#' @keywords internal
.ctGenerateJulia <- function(model, n.subjects, times, project = NULL,
  quiet = FALSE) {

  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    stop("Generating with backend='julia' needs the JuliaConnectoR package.",
      call. = FALSE)
  }
  model <- .ctGenerateResolveFree(model, quiet = quiet)
  # Individual differences are not generated here yet. Silently producing data
  # with no between-subject variation from a model that asks for it would be a
  # wrong answer rather than a missing feature, so it is refused.
  varying <- !is.null(model$pars$indvarying) && any(model$pars$indvarying)
  if (varying) {
    stop("Generating from a model with individually varying parameters is not ",
      "supported through backend='julia' yet: the specification does not carry ",
      "the population distribution generation would have to draw from. Set ",
      "indvarying to FALSE to generate a fixed-effects dataset, or use ",
      "backend='r'.", call. = FALSE)
  }

  skeleton <- .ctGenerateSkeleton(model, n.subjects, times)
  spec <- .ctJuliaPrepare(skeleton, model, project = project, priors = FALSE,
    intoverpop = "augmented")
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))

  npar <- suppressWarnings(max(c(0L, as.integer(spec$parameter_table$parnumber)),
    na.rm = TRUE))
  # A zero-length vector deadlocks the JuliaConnectoR bridge, and a fully fixed
  # model legitimately has no free parameters, so a single unread element is
  # sent instead of nothing.
  raw <- if (npar < 1L) 0 else numeric(npar)

  nmanifest <- length(model$manifestNames)
  base <- matrix(stats::rnorm(nmanifest * nrow(skeleton)), nmanifest,
    nrow(skeleton))
  # `ctsem_generate` returns the draws alongside the row likelihoods, and hands
  # them back manifest-major.
  generated <- .ctBackendGenerate(handle, raw, base)
  drawn <- if (is.list(generated)) generated$Y else generated
  drawn <- as.numeric(drawn)
  if (length(drawn) != nmanifest * nrow(skeleton)) {
    stop("The engine returned ", length(drawn), " generated values where ",
      nmanifest * nrow(skeleton), " were expected.", call. = FALSE)
  }
  values <- t(matrix(drawn, nrow = nmanifest))
  skeleton[, model$manifestNames] <- values
  # A matrix, because that is what `ctGenerate`'s own path returns and callers
  # index it positionally. Returning a data frame here would be tidier and would
  # break them.
  as.matrix(skeleton)
}
