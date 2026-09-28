# How much of each random effect a subject's own data determine.
#
# A model can be identified in principle and still ask its data for more than
# they hold. The case this is for: AnomAuth with a random CINT and a random
# DRIFT, 800 subjects with three to five waves of one indicator. Each subject's
# few observations barely constrain that subject's own rate of change, so the
# population sd of drift is barely determined, more than half the subjects'
# likelihoods go convex in their own effect, and the Laplace objective grows a
# spurious maximum 21 nats above the exact likelihood. Every one of those
# symptoms surfaced at the end of a five-minute fit, as a flat direction with
# five parameters along it, and none of them named the cause.
#
# The cause is measurable directly. For each subject and effect,
#
#     1 - posterior variance / population variance
#
# is the share of the subject's effect its own data determine: the complement
# of what population pharmacokinetics calls shrinkage. The engine computes it
# from what each objective already builds -- the unit curvature on the Laplace
# route, the filter's own covariance on the augmented one -- in one pass (see
# ctsem_effect_information in the engine). The median over subjects is the
# summary.
#
# Where it is taken: at the estimate the fit reports, after the quadrature
# correction when that moved it, beside the rest of the identifiability
# report. The plan was earlier -- at the start, so a long fit could say so
# before it spent its time -- and both earlier points were measured and
# neither works (dev/optimbench cells, local, 2026-09-27):
#
#   * At the start. The share depends on the population sd as much as on the
#     data: with per-subject information h about an effect and population sd
#     s it is s^2 h / (1 + s^2 h). At the starting values s is the start's and
#     the rest of the vector a guess, and the share then describes the start:
#     there the certified cells' random drifts read 0.04 (acnonlin) and -0.01
#     (gN1), at or below AnomAuth's 0.01.
#   * Where the optimiser stopped, before the correction. That is the Laplace
#     optimum, which can sit exactly where the approximation is least
#     trustworthy: on gD3, certified and 0.01 to 0.08 nats from its best, the
#     drift share there was -0.04 with 37 of 60 subjects convex in it, and
#     after the correction (at the best-known point) 0.47 with none.
#
# The message therefore comes at the end of the fit, with the other
# identifiability findings, and it is the cause those do not name.
#
# And why a second point when a population sd comes out small. At an sd near
# zero the share is zero by construction, whatever the data: AnomAuth's best
# point has both of its sds there, and its subjects' data pin each subject's
# CINT tightly -- no individual differences, seen clearly -- while they say
# almost nothing about each subject's drift. The share at the estimate cannot
# tell those apart; the share at a spread large enough to see can. So an
# effect whose sd is below the spread every fit starts from (raw zero, sd
# log1p_exp(-1) times its sdscale), and whose share at the estimate would
# flag it, is measured again with that sd raised to the starting spread and
# everything else at the estimate, and it is that share that decides. The
# record keeps both. An effect whose share already clears the bar at its own
# small sd is not measured again: raising the sd raises the share (for a
# linear effect exactly, s^2 h / (1 + s^2 h) with h fixed), so the second pass
# could not change the verdict and would double the check's cost on most fits,
# whose intercept sds come out small.

# The measure at one raw parameter vector, for a prepared specification.
#
# NULL when the model has no random effects or the engine cannot evaluate the
# point: this is a report, and nothing about a fit depends on it succeeding.
# `point` is the words the record and the message use for where it was taken.
#' @keywords internal
.ctEffectInformation <- function(spec, values, point) {
  started <- proc.time()[["elapsed"]]
  values <- as.numeric(values)
  at <- .ctEffectEvaluate(spec, values)
  if (is.null(at)) return(NULL)
  table <- at$table
  # Effects whose population sd is below the starting spread, whose sd is a
  # coordinate with the standard transform (raw zero is that spread there),
  # and whose share at the estimate would flag them.
  raise <- which(!is.na(table$sdindex) & table$sdindex <= length(values) &
    values[pmax(1L, table$sdindex)] < 0 & !(table$determined >=
      .ctEffectThresholds()$determined))
  table$reference <- NA_real_
  table$referencesd <- NA_real_
  if (length(raise)) {
    shifted <- values
    shifted[table$sdindex[raise]] <- 0
    again <- .ctEffectEvaluate(spec, shifted)
    if (!is.null(again) && nrow(again$table) == nrow(table)) {
      table$reference[raise] <- again$table$determined[raise]
      table$referencesd[raise] <- again$table$popsd[raise]
    }
  }
  table$weak <- .ctEffectWeak(table)
  table$sdindex <- NULL
  out <- list(point = point, values = values, route = at$route,
    table = table, groups = at$groups, thresholds = .ctEffectThresholds())
  if (!is.null(at$mins)) {
    # Per unit, as `fit$laplace$conditioning` reports it at the estimate: an
    # eigenvalue of the inner curvature below the prior's is a unit whose
    # likelihood is convex in its effects.
    out$units <- .ctJuliaLaplaceConditioning(at$mins)
    out$units$n <- length(at$mins)
  }
  out$secs <- proc.time()[["elapsed"]] - started
  out
}

# One engine pass, read into one row per effect.
#' @keywords internal
.ctEffectEvaluate <- function(spec, values) {
  laplace <- !is.null(spec$laplace)
  if (laplace) {
    levels <- spec$laplace$levels
    if (!length(levels) || !any(vapply(levels, function(l)
      isTRUE(l$nrandom > 0L), logical(1)))) return(NULL)
  } else {
    effects <- spec$random_effects
    if (is.null(effects) || !NROW(effects)) return(NULL)
    effects <- effects[effects$type %in% "sd", , drop = FALSE]
    if (!nrow(effects)) return(NULL)
  }
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- .ctJuliaModule(spec$project)
  result <- try(.ctBackendJuliaValue(module$ctsem_effect_information(
    .ctJuliaObjective(handle), .ctJuliaNumericVector(values))), silent = TRUE)
  if (inherits(result, "try-error") || is.null(result$determined)) return(NULL)
  determined <- result$determined
  if (!is.list(determined)) determined <- list(determined)
  rows <- list()
  groups <- list()
  mins <- NULL
  if (laplace) {
    mins <- as.numeric(result$min_eigenvalue)
    for (l in seq_along(levels)) {
      level <- levels[[l]]
      k <- as.integer(level$nrandom)
      if (!k) next
      d <- matrix(as.numeric(as.matrix(determined[[l]])), ncol = k)
      colnames(d) <- as.character(level$param)
      groups[[as.character(level$name)]] <- d
      popvar <- as.numeric(result$popvar[[l]])
      # What switches this level's effects off: `indvarying` for subjects,
      # `indvarying_<level>` above them.
      column <- if (l == 1L || is.null(spec$model)) "indvarying" else
        .ctJuliaLevelColumn(spec$model, l)
      unit <- if (l == 1L) "subject" else as.character(level$name)
      # A reduced-rank level has loadings rather than sds, and raw zero is not
      # a spread there but a saddle.
      sdindex <- as.integer(level$sd_index)
      if (length(sdindex) != k) sdindex <- rep(NA_integer_, k)
      for (j in seq_len(k)) rows[[length(rows) + 1L]] <- .ctEffectRow(
        level$name, level$param[j], sqrt(popvar[j]), d[, j], unit, column,
        sdindex[j])
    }
  } else {
    d <- as.matrix(determined[[1L]])
    d <- matrix(as.numeric(d), nrow = nrow(d))
    if (ncol(d) != nrow(effects)) return(NULL)
    colnames(d) <- as.character(effects$param)
    name <- .ctJuliaOr(.ctFitModelObject(list(model_spec = spec))$subjectIDname,
      "subject")
    groups[[name]] <- d
    # The filter starts each effect's state at its population variance, in
    # state units; `scale` is the factor between those and the raw parameter.
    priorsd <- sqrt(as.matrix(result$priorvar))
    for (j in seq_len(ncol(d))) {
      scale <- as.numeric(effects$scale[j])
      if (!is.finite(scale) || scale == 0) scale <- 1
      rows[[length(rows) + 1L]] <- .ctEffectRow(name, effects$param[j],
        stats::median(priorsd[, j], na.rm = TRUE) / abs(scale), d[, j],
        "subject", "indvarying",
        .ctEffectAugmentedSdIndex(spec, effects$parameter[j]))
    }
  }
  if (!length(rows)) return(NULL)
  list(table = do.call(rbind, rows), groups = groups, mins = mins,
    route = if (laplace) "laplace" else "augmented")
}

# The raw index of an augmented effect's population sd, when that coordinate
# has the sd transform (`log1p_exp(2 * param - 1)`, raw zero the starting
# spread) -- not when the sd is fixed, nor under a factor construction, whose
# diagonal is linear and whose raw zero is a saddle rather than a spread.
#' @keywords internal
.ctEffectAugmentedSdIndex <- function(spec, index) {
  index <- suppressWarnings(as.integer(index))
  if (!length(index) || is.na(index)) return(NA_integer_)
  table <- spec$parameter_table
  row <- which(!is.na(table$parnumber) & table$parnumber == index)
  transform <- if (length(row)) as.character(table$transform[row[1L]]) else ""
  if (!grepl("log1p_exp", transform, fixed = TRUE)) return(NA_integer_)
  index
}

#' @keywords internal
.ctEffectRow <- function(level, effect, popsd, determined, unit, column,
  sdindex) {
  finite <- determined[is.finite(determined)]
  data.frame(level = as.character(level), effect = as.character(effect),
    popsd = as.numeric(popsd), groups = length(finite),
    # The median over groups: what a typical subject's own data determine.
    determined = if (length(finite)) stats::median(finite) else NA_real_,
    # The share of groups whose posterior is wider than the population, by
    # more than 1%: determined below zero, a likelihood convex in the effect.
    # The margin keeps out the sign of rounding at a collapsed population sd,
    # where every share is within 1e-8 of zero and a third of AnomAuth's
    # subjects came out on the negative side of it.
    widened = if (length(finite)) mean(finite < -0.01) else NA_real_,
    unit = as.character(unit), switch = as.character(column),
    sdindex = as.integer(sdindex), stringsAsFactors = FALSE)
}

# The rule, and where it was set.
#
# An effect is reported as weakly informed when the median share its groups'
# own data determine is below `determined` -- at the starting spread when its
# population sd came out below it and the share there was needed
# (`reference`, starred below), and at the estimate otherwise. Set on the
# optimiser bench (dev/optimbench/cells.R, the baseline grid's seeds; verdicts
# from review/bench/2026-09-26-baseline.csv), from whole default fits of every
# cell with random effects, cores = 1, local, 2026-09-27:
#
#   cell (seeds)            effects: share at the estimate       flagged
#   AnomAuth S1             drift -0.00*, cint 1.00*               drift
#   AnomAuth S2             drift -0.00*, cint 1.00*               drift
#   gD1 (3)                 drift -0.00*, cint 0.98*               drift
#   gB8 (3)                 T0m 0.99*, diff 0.47, cint 1.00*       --
#   gD3 (3)                 drift 0.40 to 0.41, cint 0.98*         --
#   gN1, gN3 (nested)       drift 0.49 and 0.84, cint 0.81-0.83*,
#                           study cint 0.83*                       --
#   gA1, gA14               T0m 0.97*, drift 0.74 and 0.77,
#                           cint 0.97-0.98*                        --
#   gB2, gC2, gC8           T0m 0.96-0.99*, diffusion or
#                           measurement error 0.59-0.73, cint 0.99-1.00*  --
#   acnonlin, jflat         drift 0.47 and 0.86 (T0m, cint 0.98-0.99*)  --
#   mvmix, both routes      cint1 0.97*, cint2 0.99*               --
#   cf_gaussian, _binary,
#   _ordinal, _mixed, both  cint 0.96-0.99*                        --
#   small, bigp, panel,
#   ordinal                 cints 0.99-1.00*                       --
#
# The three flagged are the cells the check exists for (AnomAuth's
# notmaximum and saturated fits, gD1's three notstationary verdicts), and
# every one reads within 0.003 of zero; every cell that certifies reads 0.40
# or more. So 0.1, a quarter of the lowest silent share. gB8 -- random
# diffusion, seeds ending 2.4 to 3.1 nats short -- is not flagged, and no bar
# could flag it without flagging gD3 (0.40) and acnonlin (0.466) too, which
# certify with lower shares than gB8's 0.474: what goes wrong on gB8 is the
# search, not the information. At the optimiser's end point, before the
# quadrature correction, every verdict is the same but gD3's (-0.04, see
# above); at the best-known points, the same as here.
#
# Not the concavity count. A unit with an eigenvalue of its inner curvature
# below one is a unit whose likelihood is convex somewhere in its effects,
# and the nested cells gN1 and gN3, identified and certified at their best
# points, have four of their eight units there -- so `fit$laplace`'s
# `below_one` is carried in the record but decides nothing.
#' @keywords internal
.ctEffectThresholds <- function() list(determined = 0.1)

#' @keywords internal
.ctEffectWeak <- function(table) {
  reference <- table[["reference"]]
  if (is.null(reference)) reference <- rep(NA_real_, nrow(table))
  informed <- ifelse(is.finite(reference), reference, table$determined)
  is.finite(informed) & informed < .ctEffectThresholds()$determined
}

# The wording, in one place, for the message a fit prints, the note under the
# population sds in `summary()`, `ctReport()` and `print.ctIdentify()`. One
# line per weakly informed effect, a phrase in the register of
# `.ctIdentifyAdvice()`: what it means and what to do.
#' @keywords internal
.ctEffectAdvice <- function(effects) {
  table <- if (is.list(effects)) effects[["table"]] else NULL
  if (!is.data.frame(table) || !nrow(table)) return(character())
  weak <- which(table$weak %in% TRUE)
  if (!length(weak)) return(character())
  levels <- length(unique(table$level)) > 1L
  vapply(weak, function(i) {
    row <- table[i, ]
    name <- if (levels) paste0(row$effect, " (", row$level, " level)") else
      row$effect
    # Both numbers, because both are points it was evaluated at: the sd the
    # fit estimated, and the one it was measured at when that was too small
    # to measure anything. Raw scale, the scale of the effect itself, which is
    # not the transformed scale `summary()` reports the spreads on; saying so
    # is what stops a reader looking for the number in that table.
    reference <- row[["reference"]]
    if (is.null(reference)) reference <- NA_real_
    how <- if (is.finite(reference)) paste0("its raw-scale population sd was ",
      "estimated at ", signif(row$popsd, 2), ", and even at ",
      signif(row$referencesd, 2), " a typical ", row$unit, "'s own data ",
      "would determine ", .ctEffectShare(reference), " of its value") else
      paste0("at its estimated raw-scale population sd of ",
        signif(row$popsd, 2), " a typical ", row$unit, "'s own data ",
        "determine ", .ctEffectShare(row$determined), " of its value, so that ",
        "sd is poorly determined")
    paste0("Individual differences in ", name, " are barely informed: ", how,
      ". Consider ", row$switch, " = FALSE for ", row$effect,
      ", or more observations per ", row$unit, ".")
  }, character(1))
}

#' @keywords internal
.ctEffectShare <- function(x) {
  if (!is.finite(x) || x < 0.005) return("almost none")
  paste0(signif(100 * x, 2), "%")
}

# Said once, when the check is taken during a fit, and nothing when every
# effect is informed.
#' @keywords internal
.ctEffectMessage <- function(effects) {
  lines <- .ctEffectAdvice(effects)
  if (length(lines)) message(paste(lines, collapse = "\n"))
  invisible(lines)
}
