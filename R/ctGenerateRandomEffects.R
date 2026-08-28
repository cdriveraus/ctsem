# Individual differences and TI-predictor effects, in generated data.
#
# The engine already draws them. That is the whole finding this file rests on,
# and it is worth stating plainly because it changes what the code has to do.
#
# With `intoverpop = 'augmented'` an individually varying parameter *becomes a
# state*: `mm` appears as row 2 of T0MEANS, and its population standard
# deviation as `julia_popcov_2_2` in T0VAR. `ctsem_generate` already draws the
# T0 state from its own distribution, so between-subject variation comes out
# with no drawing machinery at all -- measured on a four-subject probe, subject
# means -1.77, 2.32, -2.53, 3.17 against a within-subject sd of 0.55.
#
# So there is no sampler to write. There is only a mapping to get right: from
# the population standard deviation a user states to the raw value the engine
# needs. That is what this file does.
#
# ## Which scale a user states it on
#
# The natural scale of the parameter, matching `ctGenerate()`'s own path, where
# `TRAITVAR`, `MANIFESTTRAITVAR` and `TIPREDEFFECT` are all natural-scale
# matrices applied directly to CINT and MANIFESTMEANS. Someone who writes
# `TRAITVAR = 0.3` there and `indvaryingsd = 0.3` here should get the same
# spread, and asking them to think in unconstrained units for one path and not
# the other would be gratuitous.
#
# The conversion is the derivative of the parameter's own transform at the
# population mean. For a mean parameter -- `10 * param` for T0MEANS,
# MANIFESTMEANS and CINT -- that derivative is the constant 10, so the mapping
# is exact and the population distribution is exactly normal on the natural
# scale. For a parameter with a nonlinear transform, a variance or a drift, it
# is a first-order approximation and the induced distribution is not normal.
# That is not a defect of this code but of the question: "the population SD of a
# parameter constrained to be positive" has no exact answer, and the delta
# method is the standard reading of it. Generation says which case it is in
# rather than leaving the approximation silent.

#' @keywords internal
.ctGenerateTransformEnv <- function() {
  # The transforms are ctsem's own strings and reference ctsem's own helpers, so
  # they are evaluated in an environment that has them and nothing else. Not the
  # global environment: a user variable named `param` would otherwise change how
  # a model generates.
  list2env(list(log1p_exp = log1p_exp), parent = baseenv())
}

# Evaluate a transform string at a raw value.
#
# Two forms occur: `10 * param` on the model's own rows, and `10 * param[3]` on
# the prepared parameter table. Both are handled by binding `param` to a vector
# long enough to index and reading the scalar out.
#' @keywords internal
.ctGenerateTransformAt <- function(transform, value, index = 1L) {
  env <- .ctGenerateTransformEnv()
  vector <- rep(as.numeric(value), max(1L, as.integer(index)))
  assign("param", vector, envir = env)
  result <- try(eval(parse(text = transform), envir = env), silent = TRUE)
  if (inherits(result, "try-error")) return(NA_real_)
  as.numeric(result)[1L]
}

# The slope of a transform, by central difference.
#
# Numerical rather than symbolic because the transforms are arbitrary strings a
# user may have written, and a symbolic differentiator would have to fail on
# whatever it did not recognise. A central difference on a smooth monotone
# function of one variable is accurate to far better than generation needs.
#' @keywords internal
.ctGenerateTransformSlope <- function(transform, at, index = 1L, step = 1e-5) {
  up <- .ctGenerateTransformAt(transform, at + step, index)
  down <- .ctGenerateTransformAt(transform, at - step, index)
  if (!is.finite(up) || !is.finite(down)) return(NA_real_)
  (up - down) / (2 * step)
}

# The raw value whose transform equals `target`.
#
# `uniroot` over a bracketing interval rather than an analytic inverse, for the
# same reason the slope is numerical: the transform is a string, and every
# transform ctsem generates is monotone, which is all uniroot needs.
#' @keywords internal
.ctGenerateTransformInvert <- function(transform, target, index = 1L,
  interval = c(-30, 30)) {
  if (!is.finite(target)) return(NA_real_)
  f <- function(x) .ctGenerateTransformAt(transform, x, index) - target
  lower <- f(interval[1L])
  upper <- f(interval[2L])
  if (!is.finite(lower) || !is.finite(upper) || lower * upper > 0) {
    # Outside what the transform can reach. A population SD below the
    # transform's floor or above anything it attains is a specification error
    # worth naming, not a value to clamp silently to the nearest endpoint.
    return(NA_real_)
  }
  result <- try(stats::uniroot(f, interval = interval, tol = 1e-10)$root,
    silent = TRUE)
  if (inherits(result, "try-error")) return(NA_real_)
  result
}

# The population mean's raw value for a varying parameter, which is where the
# transform's slope has to be taken.
#' @keywords internal
.ctGenerateMeanIndex <- function(spec, param) {
  table <- spec$parameter_table
  rows <- which(!is.na(table$param) & as.character(table$param) == param &
    !is.na(table$parnumber) & table$parnumber > 0)
  if (!length(rows)) return(NA_integer_)
  as.integer(table$parnumber[rows[1L]])
}

#' @keywords internal
.ctGenerateMeanRaw <- function(spec, raw, param) {
  index <- .ctGenerateMeanIndex(spec, param)
  if (is.na(index) || index < 1L || index > length(raw)) return(0)
  raw[index]
}

# Fill the raw vector's population-SD and TI-effect entries from the model.
#
# Anything the model does not state is left where it is: the raw vector arrives
# at zero, which is the centre of ctsem's prior on the population standard
# deviation, so an unstated spread is the one that prior considers typical.
# Deliberately the prior's centre rather than a draw from it -- a draw would
# make two generations from one specification differ in a way `set.seed()` hides
# rather than controls, and the point of the fallback is that an unspecified
# model still produces something sensible, not something random.
#' @keywords internal
.ctGenerateRandomRaw <- function(model, spec, raw, quiet = FALSE) {
  effects <- spec$random_effects
  pars <- model$pars
  stated <- character()
  approximate <- character()
  unreachable <- character()

  # The population mean first, because the transform's slope is taken there and
  # the answer depends on where "there" is. A varying parameter stays free
  # through preparation -- fixing it would remove the random effect entirely --
  # so its mean arrives as an intention recorded by `.ctGenerateResolveFree()`
  # rather than as a value on the row.
  intended <- attr(model, "ctGenerateMeans")
  for (param in names(intended)) {
    row <- which(!is.na(pars$param) & as.character(pars$param) == param)
    if (!length(row)) next
    index <- .ctGenerateMeanIndex(spec, param)
    if (is.na(index)) next
    transform <- as.character(pars$transform[row[1L]])
    # Inverting the *model's* transform, which is the raw-to-natural map
    # whether the augmented table applies the scaling on this row or downstream
    # of it -- the two layouts differ and this is what they agree on.
    value <- if (is.na(transform) || !nzchar(transform)) intended[[param]] else
      .ctGenerateTransformInvert(transform, intended[[param]])
    if (is.finite(value)) raw[index] <- value
  }

  if (!is.null(effects) && nrow(effects) && !is.null(pars$indvaryingsd)) {
    sds <- effects[as.character(effects$type) == "sd", , drop = FALSE]
    for (i in seq_len(nrow(sds))) {
      param <- as.character(sds$param[i])
      row <- which(!is.na(pars$param) & as.character(pars$param) == param)
      if (!length(row)) next
      target <- suppressWarnings(as.numeric(pars$indvaryingsd[row[1L]]))
      if (!is.finite(target)) next
      if (target < 0) {
        stop("indvaryingsd for '", param, "' is negative. A population ",
          "standard deviation cannot be.", call. = FALSE)
      }
      index <- as.integer(sds$parameter[i])
      if (is.na(index) || index < 1L || index > length(raw)) next

      # The requested spread is on the parameter's natural scale; the engine
      # holds it on the raw scale. The transform's slope at the population mean
      # converts between them.
      transform <- as.character(pars$transform[row[1L]])
      slope <- if (is.na(transform) || !nzchar(transform)) 1 else
        .ctGenerateTransformSlope(transform,
          .ctGenerateMeanRaw(spec, raw, param))
      if (!is.finite(slope) || slope == 0) slope <- 1
      rawsd <- target / abs(slope)

      # And the population SD is itself transformed, so the raw entry is the
      # inverse of *that* transform at the raw-scale spread.
      table <- spec$parameter_table
      cell <- which(!is.na(table$parnumber) &
        as.integer(table$parnumber) == index)
      sdtransform <- if (length(cell)) as.character(table$transform[cell[1L]]) else NA
      value <- if (is.na(sdtransform) || !nzchar(sdtransform)) rawsd else
        .ctGenerateTransformInvert(sdtransform, rawsd, index)
      if (!is.finite(value)) {
        unreachable <- c(unreachable, sprintf("%s=%s", param, format(target)))
        next
      }
      raw[index] <- value
      stated <- c(stated, sprintf("%s=%s", param, format(target)))
      # Linear means exact. Anything else is the delta method, and saying which
      # is the difference between a documented approximation and a silent one.
      if (!.ctGenerateTransformIsLinear(transform)) {
        approximate <- c(approximate, param)
      }
    }
  }

  raw <- .ctGenerateTiEffectRaw(model, spec, raw)

  if (!quiet && length(unreachable)) {
    warning("A population standard deviation was outside the range its ",
      "transform can produce and was left at its default: ",
      paste(unreachable, collapse = ", "), ".", call. = FALSE)
  }
  if (!quiet && length(approximate)) {
    message("Population spread for ", paste(approximate, collapse = ", "),
      " is set to first order: these parameters have nonlinear transforms, so ",
      "the requested standard deviation is matched at the population mean ",
      "rather than exactly, and the induced population distribution is not ",
      "normal on the natural scale. Measure what was produced if it matters.")
  }
  if (!quiet && length(stated)) {
    message("Population standard deviations set for generation: ",
      paste(stated, collapse = ", "), ".")
  }
  raw
}

# Whether a transform is linear in its parameter, tested rather than parsed: a
# straight line is the one case where a stated natural-scale spread comes back
# exactly, and it covers every mean parameter ctsem writes.
#' @keywords internal
.ctGenerateTransformIsLinear <- function(transform) {
  if (is.na(transform) || !nzchar(transform)) return(TRUE)
  at <- c(-1, 0, 1)
  slopes <- vapply(at, function(x) .ctGenerateTransformSlope(transform, x),
    numeric(1))
  all(is.finite(slopes)) &&
    max(abs(slopes - slopes[1L])) <= 1e-8 * max(1, abs(slopes[1L]))
}

# TI-predictor effects, same idea and simpler: the coefficient multiplies a
# predictor to shift a parameter, so it lives on the parameter's raw scale and
# converts by the same slope.
#' @keywords internal
.ctGenerateTiEffectRaw <- function(model, spec, raw) {
  effects <- spec$ti_effects
  if (is.null(effects) || !nrow(effects)) return(raw)
  names <- model$TIpredNames
  if (!length(names)) return(raw)
  pars <- model$pars
  # The parameter index in `ti_effects` counts the varying parameters in the
  # order `random_effects` lists their standard deviations.
  order <- spec$random_effects
  order <- if (is.null(order)) character() else
    as.character(order$param[as.character(order$type) == "sd"])
  for (i in seq_len(nrow(effects))) {
    parameter <- as.integer(effects$parameter[i])
    predictor <- as.integer(effects$predictor[i])
    index <- as.integer(effects$coefficient[i])
    if (is.na(index) || index < 1L || index > length(raw)) next
    if (is.na(parameter) || parameter < 1L || parameter > length(order)) next
    if (is.na(predictor) || predictor < 1L || predictor > length(names)) next
    param <- order[parameter]
    column <- paste0(names[predictor], "_effectsize")
    if (is.null(pars[[column]])) next
    row <- which(!is.na(pars$param) & as.character(pars$param) == param)
    if (!length(row)) next
    target <- suppressWarnings(as.numeric(pars[[column]][row[1L]]))
    if (!is.finite(target)) next
    transform <- as.character(pars$transform[row[1L]])
    slope <- if (is.na(transform) || !nzchar(transform)) 1 else
      .ctGenerateTransformSlope(transform, .ctGenerateMeanRaw(spec, raw, param))
    if (!is.finite(slope) || slope == 0) slope <- 1
    raw[index] <- target / abs(slope)
  }
  raw
}
