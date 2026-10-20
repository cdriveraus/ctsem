# Smooth effects of time-independent predictors --------------------------------
#
# `ctModel(TIpredNames = ~ s(age, k = 8) + female)`: a smooth term is a
# penalised spline of a person-level covariate, acting on each parameter it
# targets through the random effects' mean, as a linear TI predictor does.
#
# The term is a set of basis columns that the rest of ctsem treats as ordinary
# TI predictors -- named in `TIpredNames`, with `_effect` columns in `pars`,
# appended to the data before fitting -- so the engine, the transforms, the
# subject parameters and generation need nothing new. What is new is how the
# columns are built and how their coefficients are estimated:
#
# - The basis is a cubic B-spline over the covariate's values across persons
#   (one value per person, not per row), knots at their quantiles, with the
#   integrated squared second derivative as its penalty. It is centred to sum
#   to zero over persons, so the parameter's own value stays the population
#   value at the average person, and reparameterised so that the penalty is
#   the identity on all but one column: that one is the linear trend, which the
#   penalty cannot see and so stays unpenalised (`<covariate>_lin`), and the
#   others (`<covariate>_s1`, ...) are the curve's departures from a line. This
#   is the mixed-model form of a P-spline.
# - The basis is built once, from the data the model first meets, and stored on
#   the model (`$TIsmooths[[covariate]]$basis`), so a fit, its predictions and
#   generation from it evaluate the same function. Beyond the range it was
#   built on the curve continues as a straight line.
# - Each (term, target parameter) pair has its own smoothing parameter: the
#   penalised coefficients are `Normal(0, 1 / lambda)` in the objective
#   (`_ctsem_penalty`), and lambda is estimated by the generalised
#   Fellner-Schall update (Wood and Fasiolo, 2017) from the fit's Hessian:
#   `lambda <- (m - lambda tr(V_bb)) / |b|^2`, with `V` the inverse of the
#   penalised information and `m` the term's penalised dimension
#   (`.ctTIsmoothLambda()`). At its fixed point the effective degrees of
#   freedom `m - lambda tr(V_bb)` equal `lambda |b|^2`. The update is iterated
#   on the fit's quadratic model, and the fit re-optimised at the new lambdas
#   until they settle (`.ctJuliaOptimiseFit()`).
#
# A per-term lambda rather than one per covariate: the targets are on
# different raw scales, and a shared one would have to scale each by
# something -- its population sd, say -- that a target without random effects
# does not have, and that would make the penalty a function of other
# parameters. mgcv gives each linear predictor's smooth its own for the same
# reason.

# The default and smallest basis dimensions. Four is a cubic polynomial with no
# interior knots: after centring, a line and two penalised columns.
.ctTIsmoothDefaultK <- 10L
.ctTIsmoothMinK <- 4L
# Where a smoothing parameter starts: the penalised coefficients' sd 1/sqrt(10),
# in raw units of the target per unit of the basis.
.ctTIsmoothLambdaStart <- 10

# TIpredNames written as a formula: each term a plain covariate name (a linear
# predictor, as a character TIpredNames gives) or `s(covariate, k = )`. Returns
# the expanded predictor names and the smooth definitions, keyed by covariate.
.ctTIpredParseFormula <- function(formula) {
  if (length(formula) != 2L) stop("TIpredNames as a formula has no left side: ",
    "write it as ~ s(age) + female.", call. = FALSE)
  rhs <- formula[[2L]]
  terms <- list()
  collect <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+"))) {
      for (part in as.list(e)[-1L]) collect(part)
    } else terms[[length(terms) + 1L]] <<- e
  }
  collect(rhs)
  names <- character()
  smooths <- list()
  for (term in terms) {
    if (is.name(term)) {
      names <- c(names, as.character(term))
      next
    }
    if (is.call(term) && identical(term[[1L]], as.name("s"))) {
      call <- match.call(function(x, k = .ctTIsmoothDefaultK) NULL, term)
      if (!is.name(call$x)) stop("s() takes a covariate's name as its first ",
        "argument: ", deparse(term), call. = FALSE)
      covariate <- as.character(call$x)
      k <- if (is.null(call$k)) .ctTIsmoothDefaultK else eval(call$k, baseenv())
      if (!is.numeric(k) || length(k) != 1L || k != round(k) ||
          k < .ctTIsmoothMinK) stop("k in ", deparse(term), " must be a whole ",
        "number of at least ", .ctTIsmoothMinK, ".", call. = FALSE)
      k <- as.integer(k)
      if (!is.null(smooths[[covariate]])) stop("s(", covariate, ") appears ",
        "twice in TIpredNames.", call. = FALSE)
      columns <- c(paste0(covariate, "_lin"),
        paste0(covariate, "_s", seq_len(k - 2L)))
      smooths[[covariate]] <- list(covariate = covariate, k = k,
        columns = columns, linear = columns[1L], penalised = columns[-1L],
        basis = NULL)
      names <- c(names, columns)
      next
    }
    stop("TIpredNames terms must be covariate names or s(covariate, k = ): ",
      deparse(term), call. = FALSE)
  }
  if (anyDuplicated(names)) stop("A covariate appears both as a linear term ",
    "and in s(), or twice: ", paste(unique(names[duplicated(names)]),
      collapse = ", "), ".", call. = FALSE)
  clash <- intersect(names(smooths), names)
  if (length(clash)) stop(clash[1L], " is both a linear term and a smooth.",
    call. = FALSE)
  list(names = names, smooths = if (length(smooths)) smooths)
}

# Which smooth a predictor column belongs to, or NA.
.ctTIsmoothOfColumn <- function(smooths, columns) {
  out <- rep(NA_character_, length(columns))
  for (s in smooths) out[columns %in% s$columns] <- s$covariate
  out
}

# The cells a smooth reaches when a cell does not say: every free parameter
# that varies between persons at some level, under `tipredDefault = TRUE`.
# Applied to rows whose `tipreds` field was not given (`explicit` FALSE).
.ctTIsmoothDefaultEffects <- function(ctspec, smooths, explicit, tipredDefault) {
  if (!length(smooths)) return(ctspec)
  iv <- toupper(trimws(as.character(ctspec$indvarying)))
  varies <- !is.na(ctspec$indvarying) & !iv %in% c("FALSE", "F", "")
  free <- !is.na(ctspec$param) & is.na(ctspec$value) &
    !grepl("[^[:alnum:]_.]", ctspec$param)
  on <- isTRUE(tipredDefault) & varies & free & !explicit
  for (s in smooths) for (column in s$columns) {
    name <- paste0(column, "_effect")
    ctspec[[name]][!explicit] <- ifelse(on[!explicit], "TRUE", "FALSE")
  }
  ctspec
}


# Basis ------------------------------------------------------------------------

# Cubic B-splines on `knots` at `x`, with their `derivs`-th derivative.
.ctTIsmoothBSpline <- function(knots, x, derivs = 0L) {
  splines::splineDesign(knots, x, ord = 4L, derivs = rep(derivs, length(x)),
    outer.ok = TRUE)
}

# The basis definition from one value per person: knots, and the k x (k - 1)
# map from the B-splines to the centred, reparameterised columns (linear
# first). See the notes at the top of the file.
.ctTIsmoothBuild <- function(x, k, covariate) {
  x <- as.numeric(x)
  distinct <- sort(unique(x))
  if (length(distinct) < k) stop("s(", covariate, ", k = ", k, ") needs at least ",
    k, " distinct values of ", covariate, " across persons, and the data have ",
    length(distinct), ". Lower k, or use ", covariate, " as a linear predictor.",
    call. = FALSE)
  range <- range(distinct)
  ninterior <- k - 4L
  interior <- if (ninterior > 0L) unname(stats::quantile(distinct,
    seq(0, 1, length.out = ninterior + 2L)[-c(1L, ninterior + 2L)], type = 7)) else
    numeric()
  knots <- c(rep(range[1L], 4L), interior, rep(range[2L], 4L))
  # The penalty, the integral of the squared second derivative over the range:
  # a cubic spline's second derivative is linear between knots, so a two-point
  # Gauss rule on each interval is exact.
  breaks <- unique(c(range[1L], interior, range[2L]))
  nodes <- c(-1, 1) / sqrt(3)
  S <- matrix(0, k, k)
  for (i in seq_len(length(breaks) - 1L)) {
    half <- (breaks[i + 1L] - breaks[i]) / 2
    if (half <= 0) next
    at <- breaks[i] + half * (nodes + 1)
    d2 <- .ctTIsmoothBSpline(knots, at, 2L)
    S <- S + half * crossprod(d2)
  }
  B <- .ctTIsmoothBSpline(knots, x)
  # Centred over persons: the constant is the parameter's own.
  Z <- qr.Q(qr(matrix(colMeans(B), ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  Sc <- crossprod(Z, S %*% Z)
  eig <- eigen((Sc + t(Sc)) / 2, symmetric = TRUE)
  m <- k - 2L
  null <- eig$values[k - 1L]
  if (abs(null) > 1e-8 * max(abs(eig$values))) stop("Internal: the smoothing ",
    "penalty of s(", covariate, ") does not have a line as its null space.",
    call. = FALSE)
  Bc <- B %*% Z
  linear <- Bc %*% eig$vectors[, k - 1L]
  slin <- sqrt(mean(linear^2))
  sign <- if (stats::cor(linear, x) < 0) -1 else 1
  pen <- eig$vectors[, seq_len(m), drop = FALSE] %*%
    diag(1 / sqrt(eig$values[seq_len(m)]), m)
  # One common scale for the penalised columns keeps the penalty the
  # identity: the smoothest of them has unit spread over persons.
  spread <- sqrt(colMeans((Bc %*% pen)^2))
  scale <- 1 / max(spread)
  transform <- Z %*% cbind(eig$vectors[, k - 1L] * sign / slin, pen * scale)
  list(knots = knots, range = range, transform = transform, k = k)
}

# The columns at covariate values `x`: linear beyond the range the basis was
# built on, continuing the curve from its end.
.ctTIsmoothEval <- function(basis, x) {
  x <- as.numeric(x)
  clamped <- pmin(pmax(x, basis$range[1L]), basis$range[2L])
  B <- .ctTIsmoothBSpline(basis$knots, clamped)
  outside <- which(x != clamped)
  if (length(outside)) B[outside, ] <- B[outside, , drop = FALSE] +
    (x[outside] - clamped[outside]) *
      .ctTIsmoothBSpline(basis$knots, clamped[outside], 1L)
  out <- B %*% basis$transform
  out
}

# The model with every smooth's basis built (from `data` where it has none) and
# the data with the basis columns added. One value of each covariate per
# person: a smooth explains differences between persons.
.ctTIsmoothData <- function(model, data) {
  smooths <- model$TIsmooths
  if (!length(smooths)) return(list(model = model, data = data))
  id <- model$subjectIDname
  for (name in names(smooths)) {
    s <- smooths[[name]]
    covariate <- s$covariate
    if (is.null(data[[covariate]])) stop("Column ", covariate, " (the covariate ",
      "of s(", covariate, ")) is not in the data.", call. = FALSE)
    x <- suppressWarnings(as.numeric(data[[covariate]]))
    perperson <- tapply(x, data[[id]], function(v) if (all(is.na(v))) 0 else
      diff(range(v, na.rm = TRUE)))
    if (any(perperson > 1e-12 * max(1, abs(x), na.rm = TRUE))) stop(covariate,
      " varies within persons; s() takes a person-level covariate, one value ",
      "per person.", call. = FALSE)
    # One value per person, from any of their rows; built from the persons who
    # have one.
    person <- as.character(data[[id]])
    value <- tapply(x, person, function(v) if (all(is.na(v))) NA_real_ else
      v[!is.na(v)][1L])
    observed <- value[!is.na(value)]
    if (!length(observed)) stop("s(", covariate, ") has no value of ", covariate,
      " for any person.", call. = FALSE)
    if (is.null(s$basis)) {
      s$basis <- .ctTIsmoothBuild(observed, s$k, covariate)
      smooths[[name]] <- s
    }
    # A person without a value is held at the mean of those with one in their
    # group at the first level above the person, else of everyone: as a
    # missing linear predictor is held at its expected value when optimising.
    absent <- names(value)[is.na(value)]
    if (length(absent)) {
      group <- if (length(model$groupIDnames))
        tapply(as.character(data[[model$groupIDnames[1L]]]), person, `[`, 1L)
      fill <- vapply(absent, function(p) {
        peers <- if (is.null(group)) NULL else
          value[!is.na(value) & group[names(value)] == group[[p]]]
        if (length(peers)) mean(peers) else mean(observed)
      }, numeric(1))
      value[absent] <- fill
      message("s(", covariate, "): ", length(absent), " person",
        if (length(absent) > 1L) "s have" else " has", " no value of ",
        covariate, " and ", if (length(absent) > 1L) "are" else "is",
        " held at the mean of ", if (is.null(group)) "those who have one" else
        "their group, or of everyone where their group has none", ".")
    }
    x <- unname(value[person])
    columns <- .ctTIsmoothEval(s$basis, x)
    for (j in seq_along(s$columns)) data[[s$columns[j]]] <- columns[, j]
  }
  model$TIsmooths <- smooths
  list(model = model, data = data)
}


# Penalty ----------------------------------------------------------------------

# One row per (smooth, target parameter): which raw coefficients it penalises
# (`index`, a list column) and its smoothing parameter. NULL without smooths.
# `lambda` comes from the model where a fit stored one for that pair, so a
# refit starts from it.
.ctTIsmoothPenalty <- function(model, ti_effects, table) {
  smooths <- model$TIsmooths
  if (!length(smooths) || is.null(ti_effects) || !nrow(ti_effects)) return(NULL)
  predictor <- model$TIpredNames[as.integer(ti_effects$predictor)]
  cells <- table[!is.na(table$parnumber), , drop = FALSE]
  cells <- cells[!duplicated(cells$parnumber), , drop = FALSE]
  parname <- stats::setNames(as.character(cells$param), cells$parnumber)
  rows <- list()
  for (s in smooths) {
    penalised <- predictor %in% s$penalised
    for (parameter in sort(unique(as.integer(ti_effects$parameter[penalised])))) {
      take <- penalised & as.integer(ti_effects$parameter) == parameter
      label <- unname(parname[as.character(parameter)])
      start <- if (label %in% names(s$lambda)) s$lambda[[label]]
      rows[[length(rows) + 1L]] <- data.frame(term = paste0("s(", s$covariate, ")"),
        covariate = s$covariate, parameter = parameter, param = label,
        m = sum(take), lambda = if (is.null(start)) .ctTIsmoothLambdaStart else
          as.numeric(start), stringsAsFactors = FALSE)
      rows[[length(rows)]]$index <- list(as.integer(ti_effects$coefficient[take]))
      rows[[length(rows)]]$linear <- as.integer(ti_effects$coefficient[
        predictor %in% s$linear & as.integer(ti_effects$parameter) == parameter])[1L]
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$group <- seq_len(nrow(out))
  out
}

# The engine's penalty vectors from the table: one entry per coefficient.
.ctTIsmoothPenaltyVectors <- function(penalty) {
  index <- unlist(penalty$index)
  list(index = as.integer(index),
    precision = as.numeric(rep(penalty$lambda, lengths(penalty$index))))
}

# The smoothing parameters' generalised Fellner-Schall fixed point on the
# quadratic model a fit's Hessian gives at its estimate. `hessian` is the log
# posterior's at `theta`, with the penalties at `penalty$lambda`, and
# `gradient` its gradient there. Returns the new lambdas, each term's
# effective degrees of freedom, and where the quadratic model puts the optimum
# at them.
.ctTIsmoothLambda <- function(penalty, hessian, theta, gradient = NULL,
  maxiter = 200L, tol = 1e-6) {
  npar <- length(theta)
  if (is.null(gradient)) gradient <- numeric(npar)
  index <- penalty$index
  old <- rep(penalty$lambda, lengths(index))
  at <- unlist(index)
  # The information without the penalties: what the data and any other prior
  # say. The penalties' own gradient comes out of `gradient` the same way.
  data <- -(hessian + t(hessian)) / 2
  data[cbind(at, at)] <- data[cbind(at, at)] - old
  datagradient <- gradient
  datagradient[at] <- datagradient[at] + old * theta[at]
  lambda <- penalty$lambda
  edf <- rep(NA_real_, length(lambda))
  for (iter in seq_len(maxiter)) {
    precision <- rep(lambda, lengths(index))
    information <- data
    information[cbind(at, at)] <- information[cbind(at, at)] + precision
    pengradient <- datagradient
    pengradient[at] <- pengradient[at] - precision * theta[at]
    V <- tryCatch(solve(information), error = function(e) NULL)
    if (is.null(V)) break
    step <- as.numeric(V %*% pengradient)
    point <- theta + step
    updated <- lambda
    for (g in seq_along(lambda)) {
      j <- index[[g]]
      edf[g] <- length(j) - lambda[g] * sum(diag(V)[j])
      size <- sum(point[j]^2)
      updated[g] <- if (size > 0) max(edf[g], 0) / size else Inf
    }
    updated <- pmin(pmax(updated, 1e-8), 1e10)
    change <- max(abs(log(updated) - log(lambda)))
    lambda <- updated
    if (change < tol) break
  }
  precision <- rep(lambda, lengths(index))
  information <- data
  information[cbind(at, at)] <- information[cbind(at, at)] + precision
  pengradient <- datagradient
  pengradient[at] <- pengradient[at] - precision * theta[at]
  V <- solve(information)
  point <- theta + as.numeric(V %*% pengradient)
  edf <- vapply(seq_along(lambda), function(g)
    length(index[[g]]) - lambda[g] * sum(diag(V)[index[[g]]]), numeric(1))
  list(lambda = lambda, edf = edf, point = point,
    hessian = -information, iterations = iter)
}

# The smooths with each target's smoothing parameter from `penalty`, by
# parameter name, for `.ctTIsmoothPenalty()` to start a refit from.
.ctTIsmoothKeepLambda <- function(smooths, penalty) {
  for (name in names(smooths)) {
    rows <- penalty$covariate %in% name
    if (any(rows)) smooths[[name]]$lambda <- stats::setNames(
      as.list(penalty$lambda[rows]), penalty$param[rows])
  }
  smooths
}

# The smooths' table a fit reports: per term and target, its smoothing
# parameter, the penalised coefficients' sd, and the effective degrees of
# freedom of the curve's departure from a line (the line adds one more).
.ctTIsmoothReport <- function(penalty, hessian = NULL, theta = NULL) {
  if (is.null(penalty) || !nrow(penalty)) return(NULL)
  edf <- penalty$edf
  if (is.null(edf) && !is.null(hessian) && !is.null(theta)) {
    information <- -(hessian + t(hessian)) / 2
    V <- tryCatch(solve(information), error = function(e) NULL)
    if (!is.null(V)) edf <- vapply(seq_len(nrow(penalty)), function(g)
      length(penalty$index[[g]]) - penalty$lambda[g] *
        sum(diag(V)[penalty$index[[g]]]), numeric(1))
  }
  data.frame(term = penalty$term, param = penalty$param, k = penalty$m + 2L,
    lambda = penalty$lambda, sd = 1 / sqrt(penalty$lambda),
    edf = if (is.null(edf)) NA_real_ else edf + 1, stringsAsFactors = FALSE)
}

# The summary's table of smooth terms: per term and target, the effective
# degrees of freedom (the line included), the penalised coefficients' sd, and
# an approximate test of the whole term, Wood's (2013): the coefficients'
# Wald statistic on the leading `round(edf)` directions of their posterior
# covariance, against a chi-square on as many. The covariance is the fit's,
# inverse of the penalised curvature. NULL without smooths.
.ctTIsmoothSummary <- function(object, digits = 3) {
  penalty <- .ctBackendSpec(object)$penalty
  if (is.null(penalty) || !nrow(penalty)) return(NULL)
  report <- .ctTIsmoothReport(penalty)
  raw <- as.numeric(.ctFitRawEstimate(object))
  hessian <- object$uncertainty$hessian
  V <- if (!is.null(hessian)) tryCatch(solve(-(hessian + t(hessian)) / 2),
    error = function(e) NULL)
  test <- t(vapply(seq_len(nrow(penalty)), function(g) {
    j <- c(penalty$linear[g], penalty$index[[g]])
    j <- j[!is.na(j)]
    if (is.null(V) || !is.finite(report$edf[g])) return(c(NA_real_, NA_real_, NA_real_))
    eig <- eigen((V[j, j] + t(V[j, j])) / 2, symmetric = TRUE)
    r <- max(1L, min(length(j), round(report$edf[g])))
    projected <- crossprod(eig$vectors[, seq_len(r), drop = FALSE], raw[j])
    statistic <- sum(projected^2 / pmax(eig$values[seq_len(r)], .Machine$double.eps))
    c(statistic, r, stats::pchisq(statistic, r, lower.tail = FALSE))
  }, numeric(3)))
  out <- data.frame(edf = report$edf, sd = report$sd, chisq = test[, 1L],
    df = test[, 2L], p = test[, 3L])
  rownames(out) <- paste0(report$term, " on ", report$param)
  out[] <- lapply(out, function(x) round(x, digits))
  out
}

# The prediction columns of every smooth, for the reports that list linear
# predictor effects and should not list these one column at a time.
.ctTIsmoothColumns <- function(model) {
  unlist(lapply(model$TIsmooths, `[[`, "columns"), use.names = FALSE)
}

# The smooth `whichTIpreds` selects in ctTIpredEffects(): by its covariate's
# name, or by any of its columns. NULL when it selects none; an error when it
# mixes a smooth with other predictors.
.ctTIsmoothSelected <- function(model, whichTIpreds) {
  smooths <- model$TIsmooths
  if (!length(smooths)) return(NULL)
  names <- if (is.character(whichTIpreds)) whichTIpreds else
    model$TIpredNames[whichTIpreds]
  hit <- unique(c(intersect(names, names(smooths)),
    stats::na.omit(.ctTIsmoothOfColumn(smooths, names))))
  if (!length(hit)) return(NULL)
  own <- c(hit, smooths[[hit[1L]]]$columns)
  if (length(hit) > 1L || !all(names %in% own)) stop("whichTIpreds selects a ",
    "smooth term with other predictors; select s(", hit[1L], ") alone by naming ",
    hit[1L], ".", call. = FALSE)
  smooths[[hit[1L]]]
}

# Each person's own value of the plotted parameters or matrix cells, at their
# covariate value: long format, `x`, `value`, `variable` as ctPlotArrayGG()
# names its lines.
.ctTIsmoothPersons <- function(fit, smooth, labels, parmatrices) {
  spec <- .ctBackendSpec(fit)
  model <- .ctBackendModel(fit)
  first <- !duplicated(spec$data[[model$subjectIDname]])
  x <- as.numeric(spec$data[[smooth$covariate]][first])
  values <- if (isTRUE(parmatrices)) {
    extracted <- suppressMessages(ctExtract(.ctTIsmoothPointFit(fit),
      subjectMatrices = TRUE))
    vapply(labels, function(label) {
      cell <- regmatches(label, regexec("^([A-Za-z]+)\\[([0-9]+),([0-9]+)\\]$", label))[[1L]]
      if (length(cell) != 4L) return(rep(NA_real_, length(x)))
      array <- extracted[[paste0("subj_", cell[2L])]]
      if (is.null(array)) return(rep(NA_real_, length(x)))
      as.numeric(array[1L, , as.integer(cell[3L]), as.integer(cell[4L])])
    }, numeric(length(x)))
  } else {
    pars <- suppressMessages(suppressWarnings(ctSubjectPars(fit)))
    vapply(labels, function(label) if (label %in% dimnames(pars)[[3L]])
      as.numeric(pars[1L, , label]) else rep(NA_real_, length(x)),
      numeric(length(x)))
  }
  values <- matrix(values, nrow = length(x))
  out <- data.frame(x = rep(x, ncol(values)), value = as.numeric(values),
    variable = factor(rep(labels, each = length(x)), levels = labels))
  # Persons without a value of the covariate have no place on its axis.
  out[is.finite(out$value) & is.finite(out$x), , drop = FALSE]
}

# The fit at its estimate alone, for a report that wants point values.
.ctTIsmoothPointFit <- function(fit) {
  fit$estimate$rawposterior <- NULL
  fit
}

# A line per target: the smooth's effective degrees of freedom and test.
.ctTIsmoothSubtitle <- function(fit, smooth) {
  table <- .ctTIsmoothSummary(fit)
  if (is.null(table)) return(NULL)
  rows <- grepl(paste0("^s\\(", smooth$covariate, "\\) on "), rownames(table))
  if (!any(rows)) return(NULL)
  paste(sprintf("%s: edf %.1f, p %s", rownames(table)[rows], table$edf[rows],
    format.pval(table$p[rows], digits = 2)), collapse = "; ")
}
