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

#' Check the ASYMPTOTES matrix argument.
#'
#' A binary indicator's response probability is `lower + (upper - lower) F`,
#' with `F` the logistic of the linear predictor. ASYMPTOTES holds `lower` and
#' `upper` as its two columns, on the same terms as any other matrix: a number
#' fixes a cell and a name frees it. Each cell is
#' the probability itself, so a fixed ceiling with a free floor, a ceiling
#' alone, or one guessing parameter shared by name across items are all just
#' cells -- none needs a mode of its own.
#'
#' One row per binary indicator, in manifest order and named by it, so the
#' matrix holds nothing that is not an asymptote. NA means no asymptote on
#' that side -- the bound itself, 0 or 1 -- and is stored as that bound, so the
#' cells downstream are numbers and names like every other matrix's. A matrix
#' with no cell left away from its bound is no matrix at all: nothing is built,
#' fitted or reported.
#'
#' A length 2 vector is the one shorthand. It fills every binary row, which is
#' the usual model: `c('guess', NA)` is one guessing parameter for the whole
#' test.
#' @noRd
.ctCheckAsymptoteMatrix <- function(ASYMPTOTES, manifesttype, manifestNames) {
  if (is.null(ASYMPTOTES) || all(is.na(ASYMPTOTES))) return(NULL)
  binary <- manifestNames[manifesttype %in% 1L]
  nb <- length(binary)
  if (!nb) stop('ASYMPTOTES applies to binary indicators ',
    '(manifesttype 1), and this model has none', call. = FALSE)
  if (!is.matrix(ASYMPTOTES)) {
    if (length(ASYMPTOTES) != 2L) stop('ASYMPTOTES must be a matrix with one ',
      'row per binary indicator and two columns (lower, upper), or a length 2 ',
      'vector applied to every binary indicator', call. = FALSE)
    ASYMPTOTES <- matrix(ASYMPTOTES, nb, 2, byrow = TRUE)
  }
  if (!identical(dim(ASYMPTOTES), c(nb, 2L))) stop('ASYMPTOTES needs one row ',
    'per binary indicator (', nb, ': ', paste(binary, collapse = ', '),
    ') and two columns (lower, upper); got ',
    paste(dim(ASYMPTOTES), collapse = ' * '), call. = FALSE)
  out <- matrix(as.character(ASYMPTOTES), nb, 2,
    dimnames = list(binary, c('lower', 'upper')))
  out[is.na(out) | trimws(out) %in% 'NA'] <- NA
  out[is.na(out[, 1]), 1] <- '0'
  out[is.na(out[, 2]), 2] <- '1'
  .ctCheckAsymptoteValues(suppressWarnings(as.numeric(out)))
  if (!.ctAsymptoteMoved(out)) return(NULL)
  out
}

#' TRUE where a cell is free or fixed away from its bound, 0 for the lower
#' column and 1 for the upper. Takes the matrix (cells as text) or a pars
#' subset (`param`, `value`, `col`).
#' @noRd
.ctAsymptoteMoved <- function(x) {
  if (is.matrix(x)) {
    num <- suppressWarnings(matrix(as.numeric(x), nrow(x), 2))
    bound <- matrix(c(0, 1), nrow(x), 2, byrow = TRUE)
    return(any(is.na(num) | num != bound))
  }
  bound <- ifelse(x$col == 1L, 0, 1)
  !is.na(x$param) | (!is.na(x$value) & x$value != bound)
}

.ctCheckAsymptoteValues <- function(values) {
  values <- values[!is.na(values)]
  if (any(values < 0 | values > 1)) stop('a fixed asymptote is a probability, ',
    'so must lie between 0 and 1 -- or NA for none', call. = FALSE)
  invisible(NULL)
}

#' For each manifest variable, the ASYMPTOTES row it reads, or 0.
#'
#' Row k of ASYMPTOTES is the k-th binary indicator, so this is the map the
#' engine indexes with. A row whose cells all sit at their bounds maps to 0
#' and takes the plain logistic path -- the same likelihood without the
#' mixture's cost.
#'
#' Read from `pars` at fit time rather than from ctModel's argument, because a
#' cell can be freed or fixed in `pars` afterwards. And decided by the model
#' rather than by the values the engine meets: a free guessing parameter at
#' zero is an ordinary place for an optimizer to be, and is still a row with
#' asymptotes.
#' @noRd
.ctAsymptoteRows <- function(pars, manifesttype) {
  out <- integer(length(manifesttype))
  if (is.null(pars)) return(out)
  a <- pars[pars$matrix %in% 'ASYMPTOTES', , drop = FALSE]
  if (!nrow(a)) return(out)
  binary <- which(manifesttype %in% 1L)
  rows <- unique(a$row[.ctAsymptoteMoved(a)])
  rows <- rows[rows <= length(binary)]
  out[binary[rows]] <- as.integer(rows)
  out
}

#' The fit-time half of `.ctCheckAsymptoteMatrix`: `pars` can be edited after
#' ctModel, and `manifesttype` too, which would leave the rows matched to the
#' wrong items.
#' @noRd
.ctCheckAsymptotePars <- function(pars, manifesttype) {
  a <- pars[pars$matrix %in% 'ASYMPTOTES', , drop = FALSE]
  if (!nrow(a)) return(invisible(NULL))
  nb <- sum(manifesttype %in% 1L)
  if (max(a$row) != nb) stop('ASYMPTOTES has ', max(a$row), ' rows but the ',
    'model has ', nb, ' binary indicators; it needs one row per binary ',
    'indicator -- rebuild the model with ctModel()', call. = FALSE)
  .ctCheckAsymptoteValues(a$value)
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
#' @noRd
.ctThresholdMatrix <- function(ncategories, manifesttype, manifestNames) {
  n <- length(manifesttype)
  ncol <- max(ncategories - 1L)
  out <- matrix(0, n, ncol,
    dimnames = list(manifestNames, paste0('threshold', seq_len(ncol))))
  for (i in seq_len(n)) {
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

#' Binary data checked against the model that declares it.
#'
#' A Bernoulli observation is 0 or 1. Any other value is a coding error, and
#' the filter would score it anyway -- a continuous variable declared binary
#' fitted to a log likelihood with no meaning and no message.
#' @noRd
.ctDataBinary <- function(datalong, ctm) {
  binary <- which(ctm$manifesttype %in% 1)
  if (!length(binary)) return(invisible(NULL))
  for (i in binary) {
    name <- ctm$manifestNames[i]
    v <- datalong[[name]]
    if (is.null(v)) next
    v <- v[!is.na(v)]
    if (!length(v)) next
    if (any(!v %in% c(0, 1))) stop('binary variable ', name, ' has values ',
      'other than 0 and 1 (range ', min(v), ' to ', max(v), '). Code it 0/1, ',
      'or model it as manifesttype 0 if it is continuous', call. = FALSE)
    if (all(v == v[1])) warning('binary variable ', name, ' takes the single ',
      'value ', v[1], ', so its probability is not identified by the data')
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
