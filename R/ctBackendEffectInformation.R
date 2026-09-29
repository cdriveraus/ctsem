# How well the data determine each random effect's population sd.
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
# Two quantities, one the trigger and one the explanation.
#
# The explanation is per group (a subject, or a group at an outer level) and
# effect: the share of the group's effect its own data determine,
#
#     d = 1 - posterior variance / population variance,
#
# the complement of what population pharmacokinetics calls shrinkage. The
# engine computes it from what each objective already builds -- the unit
# curvature on the Laplace route, the filter's own covariance on the augmented
# one -- in one pass (ctsem_effect_information in the engine).
#
# The trigger is what all the groups together carry about the population sd,
# which a share cannot say: it grows with the number of groups, so 800
# subjects each determining little of their own effect can still pin a
# population sd down, and 20 each determining half of theirs may not. In the
# normal model of a random effect -- group g's data carry information h_g about
# its value, the population sd is s -- the share is d_g = s^2 h_g / (1 + s^2
# h_g), the group's data-based estimate of its value varies about the
# population mean with variance s^2 + 1/h_g = s^2 / d_g, and the Fisher
# information about s^2 is sum_g d_g^2 / (2 s^4). So
#
#     se(log s) = 1 / sqrt(2 n),    n = sum_g d_g^2,
#
# where `information`, n, is the number of groups whose data determined their
# own values outright (d = 1) that would carry the same information about s:
# for G such groups this is the chi-square's 1 / sqrt(2 G). A group whose
# likelihood is convex in its effect (d < 0) is counted as carrying none.
#
# It is the part of the sd's imprecision the groups' own data explain, and
# not the fit's standard error: it holds the rest of the model fixed, and it
# reads each group's likelihood through its curvature at the mode. Beside the
# fit's own Hessian on the bench (the table at `.ctEffectThresholds()`) that
# Hessian's se(log s) ran from 0.38 of it (a random measurement error, whose
# unit likelihoods are far from Gaussian) to 3.2 times it (a random diffusion
# trading off against the model's other parameters). A trade-off with another
# parameter is the identifiability report's to name, not this one's.
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
# zero every share is zero by construction, whatever the data, so the
# information is too: AnomAuth's best point has both of its sds there, and its
# subjects' data pin each subject's CINT tightly -- no individual differences,
# seen clearly -- while they say almost nothing about each subject's drift.
# The estimate cannot tell those apart; a spread large enough to see can. So
# the effects whose sd is below the spread every fit starts from (raw zero, sd
# log1p_exp(-1) times its sdscale), and whose information at the estimate
# would flag them, are measured again together with those sds raised to the
# starting spread and everything else at the estimate, and it is the
# information there that decides. The record keeps both. An effect whose
# information already clears the bar at its own small sd is not measured
# again: raising the sd raises every share (for a linear effect exactly, s^2 h
# / (1 + s^2 h) with h fixed), so the second pass could not change the verdict
# and would double the check's cost on most fits, whose intercept sds come out
# small. Together, not one at a time, because the effects share each group's
# data: raised alone, AnomAuth's drift took the subjects' differences in level
# that their CINT carries once its sd is raised as well, and read 0.79 groups
# of information against 6.5e-5 with both raised.

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
  # and whose information at the estimate would flag them.
  raise <- which(!is.na(table$sdindex) & table$sdindex <= length(values) &
    values[pmax(1L, table$sdindex)] < 0 & !(table$information >=
      .ctEffectThresholds()$information))
  table$reference <- NA_real_
  table$referenceinformation <- NA_real_
  table$referencesd <- NA_real_
  if (length(raise)) {
    shifted <- values
    shifted[table$sdindex[raise]] <- 0
    again <- .ctEffectEvaluate(spec, shifted)
    if (!is.null(again) && nrow(again$table) == nrow(table)) {
      table$reference[raise] <- again$table$determined[raise]
      table$referenceinformation[raise] <- again$table$information[raise]
      table$referencesd[raise] <- again$table$popsd[raise]
    }
  }
  table$weak <- .ctEffectWeak(table)
  table$sdindex <- NULL
  out <- list(point = point, values = values, route = at$route,
    table = table, groups = at$groups, thresholds = .ctEffectThresholds(),
    ranks = .ctEffectRanks(spec, at$route, table))
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
    name <- .ctEffectSubjectLevel(spec)
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

# The subject level's name, as `poprank` and the records call it.
#' @keywords internal
.ctEffectSubjectLevel <- function(spec) {
  .ctJuliaOr(.ctFitModelObject(list(model_spec = spec))$subjectIDname,
    "subject")
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
    # What the groups together carry about the population sd, in groups that
    # determined their own values outright; see the top of this file.
    information = sum(pmax(finite, 0)^2),
    # The share of groups whose posterior is wider than the population, by
    # more than 1%: determined below zero, a likelihood convex in the effect.
    # The margin keeps out the sign of rounding at a collapsed population sd,
    # where every share is within 1e-8 of zero and a third of AnomAuth's
    # subjects came out on the negative side of it.
    widened = if (length(finite)) mean(finite < -0.01) else NA_real_,
    unit = as.character(unit), switch = as.character(column),
    sdindex = as.integer(sdindex), stringsAsFactors = FALSE)
}

# Per level: the rank of its population covariance now, and whether a lower
# one can be asked for with `poprank`, so that the advice offers only what
# `ctFit()` would accept. Not below 1 (a rank of zero is no variation, said by
# removing the effects); not where a RAWPOPVAR cell is stated (a reduced rank
# needs the covariance free, `.ctPopRegressionRawPopVarStated()`); and on the
# augmented route not with an individually varying T0MEANS, whose carrier is
# the latent itself (see poprank in ?ctFit). A level's rank is also capped by
# its group count, which a lower rank never exceeds.
#' @keywords internal
.ctEffectRanks <- function(spec, route, table) {
  model <- spec$model
  stated <- tryCatch(length(.ctPopRegressionRawPopVarStated(model)) > 0L,
    error = function(e) TRUE)
  if (identical(route, "laplace")) {
    levels <- spec$laplace$levels
    rows <- lapply(levels, function(l) {
      k <- as.integer(l$nrandom)
      data.frame(level = as.character(l$name), effects = k,
        rank = as.integer(.ctJuliaOr(l$rank, k)), reducible = !stated,
        stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }
  k <- sum(table$level == table$level[1L])
  rank <- suppressWarnings(as.integer(model$popregression$rank))
  if (!length(rank) || is.na(rank)) rank <- k
  t0 <- any(model$pars$matrix %in% "T0MEANS" & model$pars$indvarying %in% TRUE)
  data.frame(level = as.character(table$level[1L]), effects = k, rank = rank,
    reducible = !stated && !t0, stringsAsFactors = FALSE)
}

# The rule, and where it was set.
#
# An effect is reported as weakly determined when its population sd rests on
# less `information` than two groups that determine their own values outright:
# se(log s) above 0.5, a 95% interval for the sd spanning more than a factor
# of seven. At the starting spread when its sd came out below it and the
# estimate failed the bar (`referenceinformation`, starred below), and at the
# estimate otherwise. Set on the optimiser bench (dev/optimbench/cells.R),
# from whole default fits of every cell with random effects (the baseline
# grid's seeds, cores = 1, dev2, juliaFit f4b9693b, 2026-09-29), and from the
# same check at the best-known points of starts.R, which gave the same
# verdicts. `information` [median share]; the last column is the fit's own
# Hessian's se(log s) against the bound 1 / sqrt(2 n):
#
#   cell (seeds)         groups  effect: information [share]      se: Hessian/bound
#   AnomAuth S1          800     drift 6.5e-5* [0.007*]              repaired
#                                cint 794* [1.00*]
#   AnomAuth S2          800     drift 1.8e-5* [0.001*], cint 794*   repaired
#   gD1 (3)              60      drift 0.01* [0.02*]                 repaired
#                                cint 27.6 [0.70]                    0.39/0.13
#   gN1, gN3 (3 each)    8       study cint 3.32, 2.81 [0.66, 0.59]  0.45/0.39, 0.43/0.42
#                        40      drift 11.8, 25.3; cint 16.4, 23.0   0.15-0.29/0.14-0.21
#   gC8 (3)              40      mvar 15.4-16.3 [0.06-0.10]          0.07/0.18
#   gB8 (3)              40      diff 19.1-20.4 [0.59-0.64]          0.37-0.51/0.16
#   gD3 (3)              60      drift 16.2-18.9 [0.51-0.56]         0.29/0.17
#   mvmix, both routes   30      cint2 1.5-1.6 [0.23] at its sd of
#                                0.013, 29.6* at the spread          1.0/0.57
#   small (2)            30      cint2 8.64 [0.54]                   0.54/0.24
#   acnonlin, jflat,     25-60   every effect 9.3 to 37 [0.47-1.00]  0.12-0.43/0.12-0.23
#   gA1, gA14, gB2, gC2,
#   ord4, cf_* (both
#   routes)
#   bigp, ordinal, panel 300-1000 every effect 141 to 705             0.04-0.13/0.03-0.06
#
# Flagged: AnomAuth's drift and gD1's, the cells the check exists for, at a
# hundredth of a group's information or less. Every other effect, all of them
# on cells that certify, carries 2.81 or more, the lowest the study level of
# the nested cells (eight studies, whose sd the fit's own curvature puts at
# se(log s) 0.43-0.45, itself imprecise). The gap is wide, so the bar is set at
# what it means rather than in the gap: 2, 1.4 times below that study level.
#
# The three cells the plan singled out. gD1 (binary indicators, eight waves)
# is flagged: each subject's data determine 1.6% of its own drift even at the
# starting spread, and its sd collapses. gB8 is not: its diffusion sd rests on
# half its 40 subjects, and what goes wrong there (two of three seeds end
# short) is the gated floor's basin switch, not the information. gC8 is not
# either, and it is the case the share alone gets wrong: the typical subject
# determines 6-10% of its own measurement error, below the 0.1 of the rule
# this replaces, which flagged two of its three seeds, while the 40 together
# carry 15-16 subjects' worth and the fit's own curvature puts se(log s) at
# 0.07. The Hessian is smaller than the bound there because each subject's
# likelihood is far from Gaussian in a variance effect, where the bound reads
# only its curvature at the mode; elsewhere it is larger, by up to 3.2 (gB8)
# where the sd trades off against the model's other parameters. Over every
# effect whose raw sd coordinate was above -1.5 (none of them collapsed), the
# ratio ran 0.38 to 3.2, median 1.5.
#
# Why not the Hessian's own standard error as the trigger: it is not there at
# the starting spread, where a collapsed sd has to be measured, and at a
# collapsed sd the fit's covariance is repaired along the flat direction and
# reports a standard error of 3e-7 for AnomAuth's drift. Nor the concavity
# count: a unit with an eigenvalue of its inner curvature below one is a unit
# whose likelihood is convex somewhere in its effects, and the nested cells,
# identified and certified, have four of their eight units there -- so
# `fit$laplace`'s `below_one` is carried in the record but decides nothing.
#
# One thing the rule does not see. At the spurious Laplace maxima the bench
# stores for AnomAuth (starts.R, anomS1_spurious and anomS2_spurious) the
# drift sd is 1.3-1.5 and its information 12.3 and 14.0: the Laplace
# approximation there credits the subjects with information the exact
# likelihood does not give them. The default fits of both cells end at a
# collapsed drift sd instead, where the rule fires; a fit that ended at one of
# those maxima would not be flagged.
#' @keywords internal
.ctEffectThresholds <- function() list(information = 2)

#' @keywords internal
.ctEffectWeak <- function(table) {
  reference <- table[["referenceinformation"]]
  if (is.null(reference)) reference <- rep(NA_real_, nrow(table))
  informed <- ifelse(is.finite(reference), reference, table$information)
  is.finite(informed) & informed < .ctEffectThresholds()$information
}

# The wording, in one place, for the message a fit prints, the note under the
# population sds in `summary()`, `ctReport()` and `print.ctIdentify()`. One
# line per weakly determined effect, a phrase in the register of
# `.ctIdentifyAdvice()`: the effect, the evidence, and what to do -- a lower
# `poprank` first where one can be asked for, since it keeps some of the
# effect's variation, and `indvarying = FALSE` after it.
#' @keywords internal
.ctEffectAdvice <- function(effects) {
  table <- if (is.list(effects)) effects[["table"]] else NULL
  if (!is.data.frame(table) || !nrow(table)) return(character())
  weak <- which(table$weak %in% TRUE)
  if (!length(weak)) return(character())
  column <- function(name) {
    x <- table[[name]]
    if (is.null(x)) rep(NA, nrow(table)) else x
  }
  reference <- column("reference")
  referenceinformation <- column("referenceinformation")
  referencesd <- column("referencesd")
  information <- column("information")
  levels <- unique(table$level)
  poprank <- .ctEffectPoprank(effects[["ranks"]], table, weak,
    multilevel = length(levels) > 1L)
  vapply(weak, function(i) {
    row <- table[i, ]
    name <- if (length(levels) > 1L) paste0(row$effect, " (", row$level,
      " level)") else row$effect
    units <- .ctEffectPlural(row$unit)
    # Where it was evaluated, in both numbers when the estimate was too small
    # to measure anything at: the sd the fit estimated, and the starting spread
    # it was measured at instead. Raw scale, the scale of the effect itself,
    # which is not the transformed scale `summary()` reports the spreads on;
    # saying so is what stops a reader looking for the number in that table.
    raised <- is.finite(referenceinformation[i])
    n <- if (raised) referenceinformation[i] else information[i]
    share <- if (raised) reference[i] else row$determined
    at <- if (raised) paste0("even at a raw-scale population sd of ",
      signif(referencesd[i], 2), " (estimated ", signif(row$popsd, 2), ")") else
      paste0("at its estimated raw-scale population sd of ",
        signif(row$popsd, 2))
    rests <- if (!is.finite(n) || n < 0.01) paste0("the ", row$groups, " ",
      units, if (raised) " would carry" else " carry", " almost no ",
      "information about it") else paste0("that sd ",
      if (raised) "would rest" else "rests", " on an effective ",
      signif(n, 2), " of the ", row$groups, " ", units)
    why <- paste0(", since a typical ", row$unit, "'s own data ",
      if (raised) "would determine " else "determine ", .ctEffectShare(share),
      " of its ", row$effect)
    remedy <- poprank[[as.character(row$level)]]
    remedy <- if (!is.null(remedy)) paste0("Consider ", remedy, ", else ",
      row$switch, " = FALSE for ", row$effect, ".") else
      paste0("Consider ", row$switch, " = FALSE for ", row$effect, ".")
    paste0("Individual differences in ", name, " are barely determined: ", at,
      ", ", rests, why, ". ", remedy)
  }, character(1))
}

# The `poprank` to suggest, per level with a weak effect: the level's rank less
# its weak effects, at least 1, and only when that is below the rank it has
# and the level can be reduced at all. Named per level on a multilevel model,
# which is how `ctFit()` reads it (an unnamed rank applies to every level).
#' @keywords internal
.ctEffectPoprank <- function(ranks, table, weak, multilevel) {
  if (!is.data.frame(ranks) || !nrow(ranks)) return(list())
  out <- list()
  for (level in unique(table$level[weak])) {
    r <- ranks[ranks$level == level, , drop = FALSE]
    if (nrow(r) != 1L || !isTRUE(r$reducible)) next
    lower <- max(1L, as.integer(r$rank) - sum(table$level[weak] == level))
    if (lower >= r$rank) next
    out[[level]] <- if (multilevel) paste0("poprank = c(", level, " = ",
      lower, ")") else paste0("poprank = ", lower)
  }
  out
}

#' @keywords internal
.ctEffectPlural <- function(unit) {
  if (grepl("[^aeiou]y$", unit)) sub("y$", "ies", unit) else paste0(unit, "s")
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
