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
# `ctModeltoNumeric()` sets every free parameter to zero, to stay consistent
# with ctsem's behaviour before commit 437bfbc0 and not break scripts that
# depend on it. Three matrices cannot take a literal zero on the diagonal
# because zero is singular there: DRIFT (a zero DRIFT is singular, so `fQinf()`
# cannot solve for the asymptotic covariance and generation errors outright,
# reproduced on a one-latent model with a free DRIFT, the model anyone writes
# first), and DIFFUSION and T0VAR (a zeroed diagonal on either gives a
# covariance that cannot be factorised, reproduced the same way once DRIFT
# alone was patched). All three diagonals therefore use the same near-zero
# convention `ctModel0DRIFT()` already applies when a fixed DRIFT diagonal is
# exactly zero: -1e-6 in continuous time, 1-1e-6 in discrete time for DRIFT,
# and a flat 1e-6 for DIFFUSION and T0VAR (both already variances, so there is
# no continuous/discrete distinction to make). Close enough to negligible to
# be indistinguishable in practice, but non-singular so generation can
# proceed. MANIFESTVAR stays exactly zero: it is not singular in the
# measurement update.
#
# These values move. When they do, every test whose generating model leaves a
# cell free gets different data under a `set.seed()` that reads as though it
# pinned everything -- which is how test-binary-binary-mix.R came to be
# characterising a dataset that no longer existed. A test that generates should
# state its generating model in full; the message below names every cell this
# function had to supply, so an incomplete one says so in its output.
#
# These are defaults for *simulation*, not estimates of anything, and the point
# is only that an underspecified model can still be generated from. Anything a
# user cares about they should set.
#' @keywords internal
.ctGenerateDefaults <- function(continuoustime = TRUE) {
  driftdiagonal <- if(isTRUE(continuoustime)) -1e-6 else 1 - 1e-6
  list(
    DRIFT = list(diagonal = driftdiagonal, offdiagonal = 0),
    DIFFUSION = list(diagonal = 1e-6, offdiagonal = 0),
    T0VAR = list(diagonal = 1e-6, offdiagonal = 0),
    MANIFESTVAR = list(diagonal = 0, offdiagonal = 0),
    LAMBDA = list(diagonal = 0, offdiagonal = 0),
    T0MEANS = list(diagonal = 0, offdiagonal = 0),
    MANIFESTMEANS = list(diagonal = 0, offdiagonal = 0),
    CINT = list(diagonal = 0, offdiagonal = 0),
    TDPREDEFFECT = list(diagonal = 0, offdiagonal = 0),
    PARS = list(diagonal = 0, offdiagonal = 0))
}

#' @keywords internal
.ctGenerateResolveFree <- function(model, quiet = FALSE) {
  defaults <- .ctGenerateDefaults(continuoustime = isTRUE(model$continuoustime))
  pars <- model$pars
  # A cell whose label is an expression rather than a parameter name is not a
  # free parameter waiting for a value -- it *is* the specification, and the
  # state-dependent forms this route exists to reach are written that way.
  # Filling it destroys them silently: measured, a LAMBDA cell written
  # `0.9 + 0.35 * eta1` was overwritten with the off-diagonal default of zero
  # and the generated indicator came back pure noise, with nothing in the data
  # or the message to say the loading had gone. So only a plain label is a
  # candidate, and a label that names a latent process or a predictor is a
  # reference to it rather than a parameter of its own.
  reserved <- c(model$latentNames, model$manifestNames, model$TDpredNames,
    model$TIpredNames)
  label <- !is.na(pars$param) & grepl('^[A-Za-z.][A-Za-z0-9._]*$', pars$param) &
    !pars$param %in% reserved
  free <- which(is.na(pars$value) & (is.na(pars$param) | label))
  # An individually varying parameter is filled like any other, because user
  # side generation does not draw random effects: `.ctGenerateFixedOnly()` has
  # already cleared every varying flag by the time the model reaches the
  # engine, so there is no carrier state whose mean this would have to leave
  # free. It used to be recorded on a `ctGenerateMeans` attribute and written
  # into the population mean slot later, which was the right thing to do while
  # the augmented layout was carrying the effects and is now just a slot that
  # nothing reads.
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
    # Every filled cell is named, including the ones filled with zero. Naming
    # only the non-zero ones made the common case unreadable: a model leaving
    # only MANIFESTVAR or T0MEANS free reported a count and no cells, so a
    # script -- or a test under a set.seed() -- generated from values it never
    # stated and had no way to notice when those values changed underneath it.
    filled <- c(filled, sprintf("%s[%d,%d]=%s", matrix_name,
      pars$row[i], pars$col[i], format(value)))
  }
  model$pars <- pars
  if (!quiet && length(filled)) {
    message(length(filled), " free parameter", if (length(filled) > 1) "s" else "",
      " had no value and were set for generation: ",
      paste(utils::head(filled, 8), collapse = ", "),
      if (length(filled) > 8) ", ..." else "",
      ". Set them in the model if they matter.")
  }
  model
}

# A long-format skeleton with the requested subjects and times and no
# observations: the shape generation fills in.
#
# The id and time columns are named as the *model* names them, not `id` and
# `time`. Hardcoding those two made the whole julia generation route
# unreachable for any model built with `ctModel(id = ...)` or `time = ...`:
# `.ctJuliaPrepare()` looks the columns up by `model$subjectIDname` and
# `model$timeName`, found neither, and died inside `order()` with "argument 1
# is not a vector", which names nothing about the actual mistake. Reproduced
# on a one-latent model with `id = 'subject'`, and again with `time = 'age'`.
#' @keywords internal
.ctGenerateSkeleton <- function(model, n.subjects, times) {
  idname <- model$subjectIDname
  timename <- model$timeName
  rows <- do.call(rbind, lapply(seq_len(n.subjects), function(i)
    stats::setNames(data.frame(i, times[[i]]), c(idname, timename))))
  # A grouping level above the subject needs its column present or preparation
  # cannot build the data at all, and every subject goes in one group: user
  # side generation ignores random effects entirely for now, so the grouping
  # cannot affect what is generated, and any other choice here would be a claim
  # about a design the caller did not state.
  for (nm in model$groupIDnames) rows[[nm]] <- 1L
  # Zero, not NA. The engine generates only where an observation exists --
  # `_generate_row!` is handed the observed indices and writes only those -- so
  # the skeleton's *missingness pattern* is the input and its values are not:
  # the generated draw overwrites the value and the filter's update conditions
  # on the draw, never on what was there. An all-NA skeleton therefore produces
  # an all-NA result, which is what it did before this line said zero.
  for (nm in model$manifestNames) rows[[nm]] <- 0
  # Time-dependent predictors at zero. One that is never non-zero has no effect
  # on the generated data, which is the honest default for a value the caller
  # did not supply; `ctGenerate`'s own path uses TDPREDMEANS and this could too
  # once the specification carries them.
  for (nm in model$TDpredNames) rows[[nm]] <- 0
  # Time-independent predictors are *drawn*, one value per subject, held
  # constant across that subject's rows -- which is what makes them time
  # independent. Zeros would have been consistent with the line above and wrong
  # here: a predictor column that is constant carries no information, so an
  # effect on it could not show up in the data and a fit to that data could not
  # identify one. Standard normal because the specification carries no
  # distribution for them (`ctGenerate`'s own path has TIPREDMEANS and
  # TIPREDVAR; a `ctModel` of type 'ct' has neither), and because a predictor on
  # that scale makes an effect size directly readable as "change per standard
  # deviation". Drawn on the R side so `set.seed()` governs them.
  for (nm in model$TIpredNames) {
    rows[[nm]] <- stats::rnorm(n.subjects)[rows[[idname]]]
  }
  rows
}

# Every random effect cleared, and named in a message.
#
# User side generation produces one dataset from the values a specification
# states, and nothing else: every subject gets the same parameters, and any
# individual differences the model declares are ignored. The alternative --
# drawing them -- needs a population distribution the specification only
# partly pins down, and the machinery that filled the gaps interpreted a
# RAWPOPVAR entry on the parameter's natural scale, which is not the scale the
# fit works in. So it draws nothing rather than drawing from a convention that
# disagrees with fitting.
#
# Individual differences in generated data come from `ctGenerateFromFit()`
# instead, which has a fit and therefore has the population distribution on the
# scale the fit itself used.
#
# Cleared rather than refused, because a fixed-effects draw from a model with
# random effects is a perfectly reasonable thing to want and is what the
# fixed values describe. It is the silence that would be wrong, so every
# parameter dropped is named.
#' @keywords internal
.ctGenerateFixedOnly <- function(model, quiet = FALSE) {
  columns <- .ctVaryingColumns(model)
  columns <- columns[columns %in% names(model$pars)]
  dropped <- .ctVaryingParams(model)
  for (cl in columns) model$pars[[cl]] <- FALSE
  # The population covariance goes with them: with nothing varying it describes
  # a distribution no longer in the model, and leaving it would let preparation
  # read a spread for a parameter that has no random effect.
  stated <- character()
  popvar <- model[["POPCOV"]]
  if (!is.null(popvar) && length(popvar)) {
    for (nm in rownames(popvar)) {
      if (is.finite(.ctModelPopCovValue(popvar[nm, nm]))) stated <- c(stated, nm)
    }
  }
  model[["POPCOV"]] <- NULL
  if (!is.null(model$matrices)) model$matrices$POPCOV <- NULL
  # Time-independent predictor effects go with them, and this has to be said
  # rather than left to be discovered. A TI effect shifts a *varying*
  # parameter -- the prepared spec indexes `ti_effects` against the order the
  # random effects are listed in -- so with nothing varying there is nothing
  # for an effect to shift, and `ti_effects` comes back empty. The predictor
  # column is still drawn and still varies between subjects, so the generated
  # data looks exactly like data with a predictor in it and carries no effect:
  # measured on a model stating `TI1=4.3`, the correlation between the subject
  # means and the predictor came out -0.19 over 60 subjects, which is noise.
  effects <- character()
  for (nm in model$TIpredNames) {
    column <- model$pars[[paste0(nm, "_effect")]]
    if (is.null(column)) next
    if (any(.ctTipredEffectActive(column))) effects <- c(effects, nm)
  }
  if (!quiet && (length(dropped) || length(effects))) {
    message("Generating from fixed values only. ",
      if (length(dropped)) paste0("Individual differences are ignored for ",
        paste(dropped, collapse = ", "),
        if (length(stated)) paste0(" (the population spread stated for ",
          paste(stated, collapse = ", "), " is unused)"), ". "),
      if (length(effects)) paste0("Effects of ",
        paste(effects, collapse = ", "),
        " are ignored too, since a time independent predictor acts on a ",
        "varying parameter: the predictor column is generated and carries no ",
        "effect. "),
      "Use ctGenerateFromFit() to generate with random effects.")
  }
  model
}

#' @keywords internal
.ctGenerateJulia <- function(model, n.subjects, times, project = NULL,
  quiet = FALSE, intoverstates = FALSE) {

  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    stop("Generating with backend='julia' needs the JuliaConnectoR package.",
      call. = FALSE)
  }
  model <- .ctGenerateFixedOnly(model, quiet = quiet)
  model <- .ctGenerateResolveFree(model, quiet = quiet)

  skeleton <- .ctGenerateSkeleton(model, n.subjects, times)
  # `intoverpop='augmented'` on a model with nothing varying augments nothing,
  # so this is the plain per-subject model -- the same thing `'none'` would
  # prepare, without `'none'`'s refusal of a model that declares no random
  # effects. Which it now never does, since `.ctGenerateFixedOnly()` cleared
  # them.
  spec <- .ctJuliaPrepare(skeleton, model, project = project, priors = FALSE,
    intoverpop = "augmented")
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))

  # Every source of a raw index, not just the parameter table. TI-predictor
  # coefficients and the Laplace block are numbered past the end of it, so
  # counting the table alone produced a raw vector two short and the engine
  # raised a BoundsError from deep inside `_materialize_subject_values!` --
  # which says nothing about the actual mistake. This matches how `ctFit()`
  # counts them.
  npar <- suppressWarnings(max(c(0L, as.integer(spec$parameter_table$parnumber),
    as.integer(spec$laplace$npar), as.integer(spec$ti_effects$coefficient)),
    na.rm = TRUE))
  # A zero-length vector deadlocks the JuliaConnectoR bridge, and a fully fixed
  # model legitimately has no free parameters, so a single unread element is
  # sent instead of nothing.
  raw <- if (npar < 1L) 0 else numeric(npar)

  nmanifest <- length(model$manifestNames)
  if (isTRUE(intoverstates)) {
    base <- matrix(stats::rnorm(nmanifest * nrow(skeleton)), nmanifest,
      nrow(skeleton))
    # `ctsem_generate` returns the draws alongside the row likelihoods, and
    # hands them back manifest-major.
    generated <- .ctBackendGenerate(handle, raw, base)
  } else {
    # The state-explicit route: sample the trajectory, then each observation
    # given the state at its row. The engine says how many latent innovations
    # the design needs -- one per state at each subject's first row and one per
    # diffusing state per bounded substep after it -- because it is the side
    # that knows how the intervals were split.
    #
    # Innovations first, then the observation deviates, so the two blocks are
    # drawn in a fixed order from one seed.
    nz <- .ctBackendStateDimension(handle)
    z <- stats::rnorm(nz)
    base <- matrix(stats::rnorm(nmanifest * nrow(skeleton)), nmanifest,
      nrow(skeleton))
    generated <- .ctBackendGenerateStates(handle, raw, z, base)
  }
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
