#' Thresholds for ordinal manifest variables.
#'
#' An ordinal variable with K categories is modelled by the cumulative logit
#'
#'     P(y <= k | eta) = inv_logit(tau_k - eta)
#'
#' so it needs K-1 thresholds, which have to increase. The THRESHOLDS matrix
#' holds them one variable per row, but it does not hold the thresholds
#' themselves: column 1 is tau_1 and every later column is the *gap* to the
#' previous threshold, constrained positive by its transform. The engine
#' accumulates them.
#'
#' Column 1 is fixed at zero. Shifting `mu` and every threshold together leaves
#' the cumulative logit unchanged, so one location among them is redundant, and
#' fixing `tau_1` puts that location in MANIFESTMEANS -- which is `indvarying`
#' where the thresholds are not. A person-level shift of an indicator's whole
#' category scale is then one random effect rather than one per threshold, which
#' is the reason for the choice. The free parameters are the gaps, so a K
#' category variable has K-2 of them plus its mean.
#'
#' Storing gaps rather than thresholds is what makes the ordering constraint
#' free. Handed K-1 unconstrained cells an optimiser will cross them, and a
#' crossed pair gives the category between them probability zero -- the
#' likelihood is -Inf and there is no gradient pointing back out. Writing
#' `tau_1 + exp(delta_2) + ...` into each cell's transform instead would keep
#' real thresholds in the matrix, but a transform reading several free
#' parameters is not something the julia adjoint's parameter layer supports,
#' and supporting it would cost every model a gradient per cell.
#'
#' The consequence for the user is only in how the summary reads: row m column
#' 1 is that variable's first threshold, and columns 2+ are gaps.
#'
#' @noRd
NULL

#' Check and normalise the asymptotes argument.
#'
#' How many asymptotes each binary variable estimates: 0 for the two parameter
#' logistic ctsem has always fitted, 1 for a lower asymptote (three parameter,
#' the guessing probability) and 2 for both (four parameter).
#'
#' A count per variable rather than a logical, because the three cases are
#' nested and an integer says which without a second argument. Zeroed for
#' every non-binary variable for the reason `.ctCheckNcategories` zeroes its
#' own: nothing downstream should have to remember to ignore it.
#' @noRd
.ctCheckAsymptotes <- function(asymptotes, manifesttype, manifestNames) {
  n <- length(manifesttype)
  if (is.null(asymptotes)) return(rep(0L, n))
  if (length(asymptotes) == 1L) asymptotes <- rep(asymptotes, n)
  if (length(asymptotes) != n) stop('asymptotes must have one entry per ',
    'manifest variable (', n, '), or be a single value', call. = FALSE)
  asymptotes <- as.integer(asymptotes)
  asymptotes[is.na(asymptotes)] <- 0L
  if (any(!asymptotes %in% 0:2)) stop('asymptotes must be 0 (two parameter ',
    'logistic), 1 (three parameter -- a lower asymptote) or 2 (four ',
    'parameter -- both)', call. = FALSE)
  wrong <- asymptotes > 0 & !manifesttype %in% 1L
  if (any(wrong)) stop('asymptotes apply to binary indicators (manifesttype ',
    '1) only. Check: ', paste(manifestNames[wrong], collapse = ', '),
    call. = FALSE)
  asymptotes[!manifesttype %in% 1L] <- 0L
  asymptotes
}

#' Check and normalise the ncategories argument.
#' @noRd
.ctCheckNcategories <- function(ncategories, manifesttype, manifestNames) {
  n <- length(manifesttype)
  ordinal <- manifesttype %in% 2
  if (is.null(ncategories)) {
    if (!any(ordinal)) return(rep(0L, n))
    stop('manifesttype 2 (ordinal) needs ncategories -- the number of ',
      'categories of each ordinal variable. Ordinal variable(s): ',
      paste(manifestNames[ordinal], collapse = ', '), call. = FALSE)
  }
  if (length(ncategories) == 1L) ncategories <- rep(ncategories, n)
  if (length(ncategories) != n) stop('ncategories must have one entry per ',
    'manifest variable (', n, '), or be a single value', call. = FALSE)
  ncategories <- as.integer(ncategories)
  ncategories[is.na(ncategories)] <- 0L
  # Only the ordinal entries mean anything; zero the rest so nothing downstream
  # has to remember to ignore them.
  ncategories[!ordinal] <- 0L
  bad <- ordinal & ncategories < 3L
  if (any(bad)) stop('an ordinal variable needs at least 3 categories -- use ',
    'manifesttype 1 for a binary variable. Check: ',
    paste(manifestNames[bad], collapse = ', '), call. = FALSE)
  ncategories
}

#' Build the THRESHOLDS matrix for a model with ordinal variables.
#'
#' One row per manifest variable and as many columns as the widest ordinal
#' variable needs. Variables with fewer categories leave their trailing cells
#' fixed at zero and the engine never reads them, which is what lets ordinal
#' variables with different category counts share one rectangular matrix.
#' @param asymptotes Per variable, as `.ctCheckAsymptotes` returns it. A
#' binary variable with any takes the first two columns: its lower asymptote
#' and then a gap, so that the engine can form `d = c + (1-c)g` and keep
#' `0 <= c < d <= 1` without a cell's transform reading another cell's
#' parameter. That is the same restriction the ordinal gaps exist for, and the
#' same matrix, rather than a second one shaped like it.
#' @noRd
.ctThresholdMatrix <- function(ncategories, manifesttype, manifestNames,
  asymptotes = NULL) {
  n <- length(manifesttype)
  if (is.null(asymptotes)) asymptotes <- rep(0L, n)
  ncol <- max(c(ncategories - 1L, ifelse(asymptotes > 0L, 2L, 0L)))
  out <- matrix(0, n, ncol,
    dimnames = list(manifestNames, paste0('threshold', seq_len(ncol))))
  for (i in seq_len(n)) {
    if (asymptotes[i] > 0L) {
      # Column 1 is the lower asymptote and column 2 the gap to the upper one.
      # A three parameter model fixes the gap at 1, which makes the upper
      # asymptote exactly 1 and costs the engine nothing to evaluate.
      out[i, 1] <- paste0('asymptote_', manifestNames[i])
      out[i, 2] <- if (asymptotes[i] >= 2L)
        paste0('asymptotegap_', manifestNames[i]) else 1
      next
    }
    if (!manifesttype[i] %in% 2) next
    # Column 1 stays at its initialised zero: see the file header. The location
    # lives in MANIFESTMEANS, so the ambiguous model -- both free, neither
    # identified -- cannot be built, and there is nothing left to warn about.
    #
    # Unconditional. A model that *also* fixes its manifest mean has made a
    # choice and it is a coherent one: the thresholds are then pinned to the
    # category scale and their centre is set from the latent side, by T0MEANS
    # and CINT. Nothing here second-guesses that.
    for (j in seq_len(ncategories[i] - 1L)) {
      if (j == 1L) next
      out[i, j] <- paste0('threshold_', manifestNames[i], '_', j)
    }
  }
  out
}

#' Category counts implied by the data, for checking a model against it.
#' @noRd
.ctDataCategories <- function(datalong, ctm) {
  ordinal <- which(ctm$manifesttype %in% 2)
  if (!length(ordinal)) return(invisible(NULL))
  for (i in ordinal) {
    name <- ctm$manifestNames[i]
    v <- datalong[[name]]
    v <- v[!is.na(v)]
    if (!length(v)) next
    if (any(v != round(v)) || any(v < 1)) stop('ordinal variable ', name,
      ' must be coded as consecutive integers from 1 -- found values outside ',
      'that (min ', min(v), '). Shift 0-based codes up by one.', call. = FALSE)
    k <- ctm$ncategories[i]
    if (max(v) > k) stop('ordinal variable ', name, ' has values up to ',
      max(v), ' but the model was built with ncategories = ', k, call. = FALSE)
    # Fewer observed categories than declared is legal but rarely intended: the
    # thresholds bounding an empty category are not identified by the data.
    missingcats <- setdiff(seq_len(k), unique(v))
    if (length(missingcats)) warning('ordinal variable ', name,
      ' has no observations in categor', if (length(missingcats) > 1)
        'ies ' else 'y ', paste(missingcats, collapse = ', '),
      '; the thresholds bounding those are not identified by the data')
  }
  invisible(NULL)
}

#' Count data checked against the model that declares it.
#'
#' A Poisson observation is a non-negative integer, and the two ways of getting
#' that wrong are worth separating. A negative value or a fraction is a coding
#' error and cannot be fitted at all. All-zero or near-constant counts fit
#' perfectly well but say the rate is not identified by that variable, which is
#' worth hearing before rather than after.
#' @noRd
.ctDataCounts <- function(datalong, ctm) {
  counts <- which(ctm$manifesttype %in% 3)
  if (!length(counts)) return(invisible(NULL))
  for (i in counts) {
    name <- ctm$manifestNames[i]
    v <- datalong[[name]]
    if (is.null(v)) next
    v <- v[!is.na(v)]
    if (!length(v)) next
    if (any(v < 0)) stop('count variable ', name, ' has negative values (min ',
      min(v), '); a count must be zero or more', call. = FALSE)
    if (any(v != round(v))) stop('count variable ', name, ' has non-integer ',
      'values; a count must be a whole number. Model it as manifesttype 0 if ',
      'it is a rate or a continuous measure', call. = FALSE)
    if (all(v == v[1])) warning('count variable ', name, ' takes the single ',
      'value ', v[1], ', so its rate is not identified by the data')
  }
  invisible(NULL)
}

#' Check and normalise the censoring limits.
#'
#' Limits are known constants, so the checking is about coherence rather than
#' estimability: a censored variable needs at least one finite limit, since a
#' variable censored nowhere is simply Gaussian, and the lower must be below the
#' upper. Non-censored variables carry infinities so that nothing downstream has
#' to remember to ignore them.
#' @noRd
.ctCheckCensorLimits <- function(censormin, censormax, manifesttype,
  manifestNames) {
  n <- length(manifesttype)
  censored <- manifesttype %in% 4
  expand <- function(x, default) {
    if (is.null(x)) return(rep(default, n))
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n) stop('censormin and censormax must have one entry per ',
      'manifest variable (', n, '), or be a single value', call. = FALSE)
    x <- as.numeric(x)
    x[is.na(x)] <- default
    x
  }
  lower <- expand(censormin, -Inf)
  upper <- expand(censormax, Inf)
  # Only the censored entries mean anything; the rest are neutralised so that
  # a limit left over from an edited model cannot quietly apply.
  lower[!censored] <- -Inf
  upper[!censored] <- Inf
  if (any(censored & !(lower < upper))) stop('censormin must be below ',
    'censormax. Check: ', paste(manifestNames[censored & !(lower < upper)],
      collapse = ', '), call. = FALSE)
  bad <- censored & !is.finite(lower) & !is.finite(upper)
  if (any(bad)) stop('a censored variable (manifesttype 4) needs at least one ',
    'finite limit in censormin or censormax -- censored nowhere is just a ',
    'Gaussian variable, which is manifesttype 0. Check: ',
    paste(manifestNames[bad], collapse = ', '), call. = FALSE)
  list(min = lower, max = upper)
}

#' Censored data checked against the limits the model declares.
#' @noRd
.ctDataCensored <- function(datalong, ctm) {
  censored <- which(ctm$manifesttype %in% 4)
  if (!length(censored)) return(invisible(NULL))
  for (i in censored) {
    name <- ctm$manifestNames[i]
    v <- datalong[[name]]
    if (is.null(v)) next
    v <- v[!is.na(v)]
    if (!length(v)) next
    lower <- ctm$censormin[i]
    upper <- ctm$censormax[i]
    # Beyond a limit is not censoring, it is a contradiction: the model says
    # the instrument could not record such a value.
    if (any(v < lower - 1e-8)) stop('censored variable ', name, ' has values ',
      'below its censormin of ', lower, ' (min ', min(v), '). A censored ',
      'value is recorded *at* the limit, not past it', call. = FALSE)
    if (any(v > upper + 1e-8)) stop('censored variable ', name, ' has values ',
      'above its censormax of ', upper, ' (max ', max(v), '). A censored ',
      'value is recorded *at* the limit, not past it', call. = FALSE)
    atlimit <- sum(v <= lower + 1e-8) + sum(v >= upper - 1e-8)
    if (atlimit == 0) warning('censored variable ', name, ' has no ',
      'observations at either limit, so the censoring never applies and the ',
      'fit is the same as manifesttype 0')
    if (atlimit == length(v)) warning('censored variable ', name, ' has every ',
      'observation at a limit, so it carries no information about the ',
      'location of the latent process beyond which side it fell')
  }
  invisible(NULL)
}

#' Ordinal thresholds, cumulated and put on the latent scale.
#'
#' The free parameters are gaps (see the file header), and a gap is not a
#' quantity anyone wants to read: `threshold_y1_2` is the distance from the
#' first threshold to the second, the first is fixed at zero, and the location
#' of the whole set is in MANIFESTMEANS. So the numbers in `popmeans` cannot be
#' compared to the latent process, to each other across items, or to anything a
#' reader knows. Three steps recover the thresholds themselves:
#'
#'   1. cumulate the row, which is what the engine does
#'      (`_ordinal_thresholds!` in binary_measurement.jl);
#'   2. subtract MANIFESTMEANS, because the link is
#'      `inv_logit(tau_k - (LAMBDA eta + MANIFESTMEANS))`, so the location
#'      lives on the manifest side;
#'   3. divide by the loading, to get from the linear predictor's units to the
#'      latent process's own.
#'
#' Done per draw rather than on the summarised gaps, so the intervals are the
#' intervals of the cumulated quantity rather than a delta method on it.
#'
#' Step 3 needs the item to load on exactly one latent. Where it loads on
#' several there is no single scale to express a threshold in, and those rows
#' are left in linear predictor units.
#'
#' The thresholds themselves are all this reports. A column expressing them in
#' the process's own spread was tried and removed: dividing by the stationary
#' sd is a z-score only for a process centred at zero, and the obvious repair
#' -- centring on `asymCINT` -- puts `-DRIFT^-1 CINT` into a number read at a
#' glance, which is unbounded as the process slows, absent for a
#' non-stationary system, and a linearisation wherever DRIFT depends on the
#' state. Neither version is a diagnostic, so there is none.
#' @noRd
.ctThresholdSummary <- function(object, flat, layout,
  digits = 3, chains = NULL) {
  model <- .ctFitModelObject(object)
  ordinal <- which(model$manifesttype %in% 2L)
  if (!length(ordinal)) return(NULL)

  arrays <- .ctBackendPopArraysFromFlat(flat, layout, .ctBackendSpec(object))
  gaps <- arrays$pop_THRESHOLDS
  location <- arrays$pop_MANIFESTMEANS
  loadings <- arrays$pop_LAMBDA
  if (is.null(gaps) || is.null(location) || is.null(loadings)) return(NULL)

  ndraws <- dim(gaps)[1L]
  values <- list()
  for (i in ordinal) {
    k <- min(max(model$ncategories[i] - 1L, 0L), dim(gaps)[3L])
    if (k < 1L) next
    g <- matrix(gaps[, i, seq_len(k)], nrow = ndraws)
    # Cumulated by a loop rather than apply(): with one free gap apply() hands
    # back a vector and the transpose that fixes the general case breaks this
    # one, silently, in the shape of the result.
    tau <- g
    if (k > 1L) for (j in 2:k) tau[, j] <- tau[, j - 1L] + g[, j]
    tau <- tau - location[, i, 1L]

    lambda <- matrix(loadings[, i, ], nrow = ndraws)
    carried <- which(colMeans(abs(lambda)) > 1e-8)
    if (length(carried) == 1L) tau <- tau / lambda[, carried]

    for (j in seq_len(k)) {
      values[[paste0(model$manifestNames[i], " ", j, "|", j + 1L)]] <- tau[, j]
    }
  }
  if (!length(values)) return(NULL)

  .ctBackendSampleSummary(do.call(cbind, values), digits = digits,
    chains = chains)
}
