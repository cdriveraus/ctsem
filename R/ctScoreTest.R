# Score tests of relaxations of a fitted model.
#
# Every relaxation tested here is an extra coordinate whose value zero is the
# fitted model: a random effect on a parameter that does not yet vary, a time
# independent covariate's effect on a parameter, a term making the drift depend
# on the latent state. The test needs only the score in that coordinate and the
# curvature blocks that involve it, evaluated at the fit -- no refit -- and
# projects the fitted parameters out of it (Neyman's C(alpha) form), so the
# test holds at an estimate that is near, not exactly at, the maximum (priors,
# a finish that stopped on gain). The design, the literature and the checks
# behind it are in CT-SEM/review/SCORETESTS-misfit-design-2026-10-09.md.
#
# Two sources of the pieces:
#
#   * Per-subject: each subject's score s_i (the engine's own per-subject
#     gradient, `ctsem_subject_gradients`, which is the marginal score on every
#     route) and its Hessian H_i, by central differences of those exact
#     gradients. A covariate's effect adds to a parameter's raw value, so its
#     score is z_i s_ij exactly; a random effect's variance has score
#     (s_ij^2 + H_i,jj) / 2 exactly (Chesher, 1984). Any number of covariates
#     is then a few matrix products in R.
#   * An extended model: the fitted model with the new terms written into it at
#     zero, prepared the way the fit was (`.ctFitReplayArgs`), its gradient and
#     its Hessian columns along the new coordinates (`ctsem_hessian_columns`).
#     State dependence goes this way, because its terms enter the filter's
#     linearisation and no per-subject identity carries them.
#
# State dependence is written into CINT rather than into DRIFT cells. The
# engine's Jacobian is derived from DRIFT * state + CINT (`ctJacobian`), so a
# term d * eta_k * eta_m in CINT[i] is the same model as DRIFT[i,k] + d * eta_m,
# and leaves every fitted parameter where and as it was. Using the latent state
# inside the filter is what keeps the test valid: the filter evaluates the term
# at the state predicted from earlier rows, and a weight built from earlier
# rows keeps each row's score mean zero, where a smoothed or updated state does
# not (scripts/scoretests-2026-10/predictable_weights.R: a true model rejected
# at 100 per cent with those weights).

#' Score tests for relaxing a fitted model
#'
#' Tests, from one fitted model and without refitting it, whether the model
#' would fit better with a relaxation it does not contain: parameters that
#' differ between subjects (an unmodelled random effect), parameters that
#' depend on a time independent covariate the model does not include, or
#' dynamics that depend on the latent state (missing nonlinearity). Each test
#' is a score (Lagrange multiplier) test: the gradient of the log likelihood
#' in the direction of the relaxation, at the fitted model, set against its
#' sampling variance with the fitted parameters projected out. Each also gives
#' a one-step estimate of the relaxation's size -- roughly what a refit
#' including it would estimate -- for screening; refit with the relaxations
#' that matter before interpreting their sizes.
#'
#' @section Heterogeneity: For each parameter that is free, does not already
#'   vary over subjects, and is not an initial condition (T0MEANS and T0VAR,
#'   whose variation T0VAR already carries), the score test that its variance
#'   over subjects is zero. A variance cannot be negative, so the test is one
#'   sided: \code{z} large and positive says the subjects' data disagree with
#'   one shared value more than chance allows. The statistic is White's
#'   information matrix test in that parameter, which is the score test for
#'   neglected heterogeneity (Chesher, 1984): it responds to any failure of
#'   the information equality in that coordinate, of which subject
#'   differences are the usual cause. \code{sd} is the one-step estimate of
#'   the subjects' standard deviation, on the raw (unconstrained) scale and,
#'   for model parameters, approximately on the parameter's own scale
#'   (\code{sd_natural}). To act on a result, set \code{indvarying = TRUE} for
#'   the parameter and refit.
#'
#' @section Covariates: For each covariate, the effect it would have on every
#'   tested parameter if added as a time independent predictor (on the raw
#'   scale, as \code{TIPREDEFFECT} acts), including on the population standard
#'   deviations of parameters that vary (a covariate moderating
#'   heterogeneity). \code{joint} tests all of a covariate's effects together
#'   and is the screen: it is valid whatever the covariate does to other
#'   parameters. Per parameter, \code{p} tests that parameter's effect given
#'   the covariate's effects on all the others (from the joint one-step
#'   estimate), and \code{p_alone} tests it as though no other parameter were
#'   affected. The two differ when a covariate that shifts one parameter --
#'   the level of a random effect, say -- shows up in another's score through
#'   their correlation; trust \code{p}. \code{estimate} is the one-step
#'   effect per unit of the covariate, raw and (\code{estimate_natural})
#'   approximately on the parameter's own scale. A covariate is entered
#'   linearly (centred), as treatment indicators for a factor, or as a
#'   natural spline with \code{df} degrees of freedom; subjects missing it
#'   contribute at the mean. To act on a result, add the covariate to the data
#'   and to \code{TIpredNames} and refit.
#'
#' @section State dependence: For each latent equation and each pair of
#'   latent states \code{k <= m}, the score test of a term \code{d * eta_k *
#'   eta_m} in that equation's drift: the second-order (Taylor) expansion of a
#'   drift that depends on the state (Luukkonen, Saikkonen & Terasvirta, 1988,
#'   in continuous time). A term \code{eta_k * eta_m} in the equation of
#'   \code{eta_i} is \code{DRIFT[i,k]} varying with \code{eta_m} -- and,
#'   when \code{k != m}, equally \code{DRIFT[i,m]} varying with
#'   \code{eta_k}; the data cannot tell those apart at this order.
#'   \code{estimate} is the one-step change in that drift coefficient per
#'   unit of the state. The terms are evaluated by the model's own filter at
#'   the state predicted from earlier observations, which keeps the test
#'   valid; a state estimate that has seen the observation being scored (the
#'   updated or smoothed state) would not. To act on a result, write the
#'   dependence into the DRIFT cell (\code{'dr11 + d * eta2'}, with \code{d}
#'   and \code{dr11} declared in \code{PARS}) and refit.
#'
#' @section Reading the results: The tests are screens, made at the fitted
#'   model: a badly misspecified model lights up many of them at once, and
#'   their one-step sizes are local approximations -- especially for drift
#'   parameters, which act through the matrix exponential. Within each family
#'   \code{p_adj} is the Benjamini-Hochberg adjustment over the tests shown.
#'   The tests use the fit's own objective, so they apply to the route it was
#'   fitted on (\code{intoverpop}); they are made at an optimum and refuse a
#'   sampled fit. Heterogeneity and covariate tests need subjects as the
#'   independent units, so they are not available for a model with a level
#'   above the subject.
#'
#' @param fit A fit from \code{\link{ctFit}} with \code{backend = 'julia'},
#'   optimised.
#' @param covariates Optional data frame of subject-level covariates to test:
#'   a column named as the fit's subject id (or \code{id}) and one column per
#'   covariate. Long data with repeated rows per subject is accepted when each
#'   covariate is constant within subject. Without an id column, one row per
#'   subject in the fit's subject order.
#' @param tests Which families to run: any of \code{'heterogeneity'},
#'   \code{'covariates'}, \code{'statedependence'}. Covariates are tested only
#'   when \code{covariates} is given.
#' @param form How a covariate enters: \code{'auto'} (indicators for a factor,
#'   character or logical covariate; linear otherwise), \code{'linear'},
#'   \code{'spline'} or \code{'factor'}. One value, or one per covariate.
#' @param df Degrees of freedom of the spline basis.
#' @param parameters Optional parameter names to restrict the heterogeneity and
#'   covariate tests to.
#' @param order Highest degree of the state-dependence terms: 2 tests each
#'   drift coefficient varying linearly with each state; 3 adds variation with
#'   products and squares of states, which is what a symmetric dependence --
#'   regulation that strengthens at both extremes, say -- shows up as, and
#'   which degree 2 cannot see.
#' @param cores Engine cores for the evaluations.
#' @param verbose Report progress.
#' @return An object of class \code{ctScoreTest}: a list with
#'   \code{heterogeneity} (a data frame, one row per parameter),
#'   \code{covariates} (one element per covariate, each with its \code{joint}
#'   test and a data frame \code{parameters}), \code{statedependence} (a data
#'   frame of terms, with \code{joint} tests per equation and overall), and
#'   \code{info} (the route, the number of subjects, the evaluation counts).
#' @references Chesher, A. (1984). Testing for neglected heterogeneity.
#'   Econometrica, 52, 865-872.
#'
#'   Luukkonen, R., Saikkonen, P., & Terasvirta, T. (1988). Testing linearity
#'   against smooth transition autoregressive models. Biometrika, 75, 491-499.
#'
#'   Sorbom, D. (1989). Model modification. Psychometrika, 54, 371-384.
#' @examples
#' \dontrun{
#' fit <- ctFit(dat, model, backend = 'julia')
#' st <- ctScoreTest(fit, covariates = subjectdata)
#' st
#' st$heterogeneity
#' }
#' @export
ctScoreTest <- function(fit, covariates = NULL,
  tests = c('heterogeneity', 'covariates', 'statedependence'),
  form = c('auto', 'linear', 'spline', 'factor'), df = 3, parameters = NULL,
  order = 2L, cores = 2L, verbose = TRUE) {
  if (!inherits(fit, 'ctJuliaFit')) stop("ctScoreTest() needs a fit from ",
    "ctFit(backend = 'julia'). Refit on that backend to use it.", call. = FALSE)
  if (!is.null(fit$sample)) stop("Score tests are made at an optimum, and this ",
    "fit was sampled. Fit with optimize = TRUE to test it.", call. = FALSE)
  if (isFALSE(fit$args$resolved$intoverstates)) stop("intoverstates = FALSE ",
    "targets the joint density of parameters and states, whose maximum is not ",
    "an estimate of anything (see ?ctFit); there is no fitted model to test.",
    call. = FALSE)
  testsgiven <- !missing(tests)
  tests <- match.arg(tests, several.ok = TRUE)
  if ('covariates' %in% tests && is.null(covariates)) {
    if (testsgiven) stop("tests includes 'covariates' but no covariates were ",
      "given.", call. = FALSE)
    tests <- setdiff(tests, 'covariates')
  }
  if (!is.null(covariates) && !'covariates' %in% tests && testsgiven)
    message("covariates given but not tested: tests does not include 'covariates'.")
  # One form for every covariate unless several are given; match.arg() with
  # several.ok would return the whole default vector for a missing argument.
  form <- if (missing(form)) 'auto' else match.arg(form, several.ok = TRUE)
  say <- function(...) if (isTRUE(verbose)) message(...)

  spec <- fit$model_spec
  est <- .ctFitRawEstimate(fit)
  ptab <- .ctScoreParameterTable(fit)
  out <- list(heterogeneity = NULL, covariates = NULL, statedependence = NULL,
    info = list(route = .ctBackendIntOverPop(spec),
      nsubjects = .ctFitNsubjects(fit), npar = length(est)))

  persubject <- intersect(tests, c('heterogeneity', 'covariates'))
  pieces <- NULL
  if (length(persubject)) {
    het <- .ctScoreSelect(ptab, parameters, 'heterogeneity')
    cov <- .ctScoreSelect(ptab, parameters, 'covariates')
    if (!'heterogeneity' %in% tests) het <- integer()
    if (!'covariates' %in% tests) cov <- integer()
    say("Score contributions: ", 2L * length(est) + 1L, " per-subject gradient passes")
    pieces <- .ctScorePieces(fit, est, rows = union(het, cov), second = het,
      cores = cores)
    if (nrow(pieces$S) != .ctFitNsubjects(fit)) stop("This model has a level ",
      "above the subject, so its independent units are not subjects. ",
      "Heterogeneity and covariate tests need subjects; use ",
      "tests = 'statedependence'.", call. = FALSE)
    if (length(het)) out$heterogeneity <- .ctScoreHeterogeneity(pieces, ptab, het)
    if (length(cov)) {
      values <- .ctScoreCovariateValues(fit, covariates)
      forms <- rep_len(form, length(values))
      out$covariates <- Map(function(z, f, name) .ctScoreCovariate(pieces, ptab,
        cov, .ctScoreBasis(z, f, df, name)), values, forms, names(values))
    }
  }
  if ('statedependence' %in% tests) {
    say("State dependence: preparing the model with the drift terms at zero")
    information <- if (!is.null(pieces)) -pieces$Hsum else NULL
    out$statedependence <- .ctScoreStateDependence(fit, est, information,
      order = order, cores = cores)
  }
  structure(out, class = 'ctScoreTest')
}

# One row per raw parameter: its name, what kind it is, the cell it lives in,
# whether it already varies over subjects, and the slope of its transform at
# the estimate (to put one-step sizes on the parameter's own scale).
.ctScoreParameterTable <- function(fit) {
  est <- .ctFitRawEstimate(fit)
  name <- as.character(.ctFitRawParNames(fit))
  kind <- ifelse(grepl('^popsd_', name), 'popsd',
    ifelse(grepl('^(rawcor_|poploading_|popunique_|rawtipredeffect_|tipredvalue_|raw\\[)', name),
      'other', 'model'))
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) base <- .ctFitModelObject(fit)
  pars <- base$pars
  at <- match(name, pars$param)
  ns <- asNamespace('ctsem')
  slope <- vapply(seq_along(name), function(k) {
    if (kind[k] != 'model' || is.na(at[k])) return(NA_real_)
    tf <- pars$transform[at[k]]
    if (is.null(tf) || is.na(tf) || !nzchar(tf)) return(1)
    f <- function(x) eval(parse(text = tf), list(param = x), ns)
    h <- 1e-6 * max(1, abs(est[k]))
    out <- try((f(est[k] + h) - f(est[k] - h)) / (2 * h), silent = TRUE)
    if (inherits(out, 'try-error') || !is.finite(out)) NA_real_ else out
  }, numeric(1))
  data.frame(index = seq_along(name), name = name, kind = kind,
    matrix = ifelse(is.na(at), NA_character_, as.character(pars$matrix[at])),
    row = ifelse(is.na(at), NA_integer_, as.integer(pars$row[at])),
    col = ifelse(is.na(at), NA_integer_, as.integer(pars$col[at])),
    indvarying = !is.na(at) & as.logical(pars$indvarying[at]) %in% TRUE,
    slope = slope, stringsAsFactors = FALSE)
}

# The raw positions a family tests: model parameters for both, minus those
# already varying and the initial conditions for heterogeneity; plus the
# population scales for covariates. `parameters` narrows either.
.ctScoreSelect <- function(ptab, parameters, family) {
  keep <- if (family == 'heterogeneity') {
    ptab$kind == 'model' & !ptab$indvarying &
      !ptab$matrix %in% c('T0MEANS', 'T0VAR')
  } else ptab$kind %in% c('model', 'popsd')
  if (!is.null(parameters)) {
    unknown <- setdiff(parameters, ptab$name)
    if (length(unknown)) stop("Unknown parameters: ", paste(unknown, collapse = ', '),
      ". The fit's are: ", paste(ptab$name, collapse = ', '), call. = FALSE)
    keep <- keep & ptab$name %in% parameters
  }
  ptab$index[keep]
}

# Each subject's score at the estimate, and its Hessian rows for the tested
# parameters, by central differences of the engine's exact per-subject
# gradients (2p + 1 passes). The same passes give, for the heterogeneity
# parameters, the second derivative of every score component along that
# parameter, which is the third derivative the heterogeneity test's projection
# needs: d H_i[j,j] / d theta_k = d^2 s_ik / d theta_j^2. A step of 1e-3 on the
# raw scale balances truncation (h^2) against the rounding of the second
# difference (eps / h^2).
.ctScorePieces <- function(fit, est, rows, second, step = 1e-3, cores = 2L) {
  p <- length(est)
  if (!length(rows)) stop("No parameters to test.", call. = FALSE)
  .ctBackendWithMaxChunks(cores, {
    S <- .ctBackendScoreMatrix(fit, est)
    n <- nrow(S)
    H <- array(0, c(n, length(rows), p))
    D2 <- array(0, c(n, p, length(second)))
    Hsum <- matrix(0, p, p)
    for (k in seq_len(p)) {
      h <- step * max(1, abs(est[k]))
      e <- est; e[k] <- est[k] + h
      plus <- .ctBackendScoreMatrix(fit, e)
      e[k] <- est[k] - h
      minus <- .ctBackendScoreMatrix(fit, e)
      column <- (plus - minus) / (2 * h)
      Hsum[, k] <- colSums(column)
      H[, , k] <- column[, rows, drop = FALSE]
      j <- match(k, second)
      if (!is.na(j)) D2[, , j] <- (plus - 2 * S + minus) / h^2
    }
  })
  Hsum <- (Hsum + t(Hsum)) / 2
  list(S = S, H = H, rows = rows, D2 = D2, second = second, Hsum = Hsum)
}

# A symmetric pseudo-inverse and its rank: a direction can carry no
# information here (two parameters whose scores coincide) without that being
# an error.
.ctScorePinv <- function(M, rtol = 1e-10) {
  M <- (M + t(M)) / 2
  e <- eigen(M, symmetric = TRUE)
  keep <- e$values > rtol * max(abs(e$values), .Machine$double.eps)
  inv <- e$vectors[, keep, drop = FALSE] %*%
    (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
  list(inv = inv, rank = sum(keep))
}

.ctScoreBH <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- stats::p.adjust(p[ok], 'BH')
  out
}

# The variance score for each parameter j, (s_ij^2 + H_i,jj) / 2 summed over
# subjects, projected on the fitted parameters. Its variance per subject is
# H_i,jj^2 / 2 -- exact when s_ij is Gaussian, as it is for a parameter that
# enters a linear Gaussian model's mean, and the score's large-sample
# behaviour otherwise -- less what estimating theta removes,
# D' I^-1 D with D the derivative of the statistic in theta.
.ctScoreHeterogeneity <- function(pieces, ptab, het) {
  S <- pieces$S
  I <- -pieces$Hsum
  Iinv <- .ctScorePinv(I)$inv
  Stot <- colSums(S)
  rows <- lapply(seq_along(het), function(a) {
    j <- het[a]
    r <- match(j, pieces$rows)
    s <- S[, j]
    Hj <- pieces$H[, r, , drop = TRUE]
    if (is.null(dim(Hj))) Hj <- matrix(Hj, nrow = nrow(S))
    hjj <- Hj[, j]
    stat <- sum(0.5 * (s^2 + hjj))
    # D = dT/dtheta. T(theta) + D' I^-1 S(theta) is T at the true theta to
    # first order wherever it is evaluated -- the C(alpha) form. Plus, because
    # D is the derivative and I the negative Hessian; the covariate test's
    # `- Idt I^-1 S` is the same with Idt = -dU/dtheta.
    D <- colSums(s * Hj) + 0.5 * colSums(pieces$D2[, , a, drop = TRUE])
    statstar <- stat + sum(D * (Iinv %*% Stot))
    V <- sum(0.5 * hjj^2) - as.numeric(t(D) %*% Iinv %*% D)
    z <- if (V > 0) statstar / sqrt(V) else NA_real_
    sdraw <- if (V > 0) sqrt(max(0, statstar / V)) else NA_real_
    data.frame(parameter = ptab$name[j], matrix = ptab$matrix[j],
      row = ptab$row[j], col = ptab$col[j], z = z,
      p = stats::pnorm(z, lower.tail = FALSE), sd = sdraw,
      sd_natural = abs(ptab$slope[j]) * sdraw, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adj <- .ctScoreBH(out$p)
  out[order(out$p), , drop = FALSE]
}

# Covariate values per subject, in the fit's subject order, as a named list.
.ctScoreCovariateValues <- function(fit, covariates) {
  covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  idname <- .ctFitModelObject(fit)$subjectIDname
  ids <- .ctFitIdMap(fit)$original
  idcol <- intersect(c(idname, 'id'), names(covariates))[1]
  vars <- setdiff(names(covariates), idcol)
  if (!length(vars)) stop("covariates has no columns besides the id.", call. = FALSE)
  inmodel <- intersect(vars, .ctFitModelObject(fit)$TIpredNames)
  if (length(inmodel)) stop("Already in the model as time independent ",
    "predictors: ", paste(inmodel, collapse = ', '), ". Test covariates the ",
    "model does not contain.", call. = FALSE)
  if (is.na(idcol)) {
    if (nrow(covariates) != length(ids)) stop("covariates needs a '", idname,
      "' column, or exactly one row per subject in the fit's subject order (",
      length(ids), ").", call. = FALSE)
    return(as.list(covariates[vars]))
  }
  key <- as.character(covariates[[idcol]])
  out <- lapply(vars, function(v) {
    x <- covariates[[v]]
    first <- match(as.character(ids), key)
    if (anyDuplicated(key)) {
      varies <- tapply(seq_along(x), key, function(r) length(unique(x[r][!is.na(x[r])])) > 1L)
      if (any(varies)) stop("Covariate '", v, "' varies within subject. Only ",
        "time independent covariates are tested here.", call. = FALSE)
      firstvalue <- tapply(seq_along(x), key, function(r) {
        r <- r[!is.na(x[r])]
        if (length(r)) r[1] else NA_integer_
      })
      first <- as.integer(firstvalue[as.character(ids)])
    }
    x[first]
  })
  names(out) <- vars
  missingids <- sum(!as.character(ids) %in% key)
  if (missingids) message(missingids, " subjects have no row in covariates and ",
    "contribute at the covariates' means.")
  out
}

# The basis a covariate enters with, centred, missing values at the mean
# (zero after centring) -- any fixed function of the covariates is a valid
# weight, so a subject without a value contributes nothing and biases nothing.
.ctScoreBasis <- function(z, form, df, name) {
  isfactor <- is.factor(z) || is.character(z) || is.logical(z)
  if (form == 'auto') form <- if (isfactor) 'factor' else 'linear'
  if (form == 'factor') {
    f <- droplevels(as.factor(z))
    if (nlevels(f) < 2L) stop("Covariate '", name, "' has fewer than two levels.",
      call. = FALSE)
    X <- matrix(0, length(z), nlevels(f) - 1L)
    for (l in seq_len(ncol(X))) X[, l] <- as.numeric(f == levels(f)[l + 1L])
    colnames(X) <- levels(f)[-1L]
  } else {
    if (isfactor) stop("Covariate '", name, "' is not numeric; use form = ",
      "'factor'.", call. = FALSE)
    z <- as.numeric(z)
    ok <- !is.na(z)
    if (length(unique(z[ok])) < 2L) stop("Covariate '", name, "' does not vary.",
      call. = FALSE)
    if (form == 'linear') {
      X <- matrix(z, ncol = 1L, dimnames = list(NULL, 'linear'))
    } else {
      X <- matrix(NA_real_, length(z), df)
      X[ok, ] <- splines::ns(z[ok], df = df)
      colnames(X) <- paste0('ns', seq_len(df))
    }
  }
  ok <- stats::complete.cases(X)
  X[ok, ] <- sweep(X[ok, , drop = FALSE], 2, colMeans(X[ok, , drop = FALSE]))
  X[!ok, ] <- 0
  attr(X, 'name') <- name
  attr(X, 'form') <- form
  attr(X, 'nmissing') <- sum(!ok)
  X
}

# The C(alpha) test of a covariate's effects on the parameters `cov`, with the
# basis G: U = G' S over those columns, and the information blocks from the
# per-subject Hessians weighted by G.
.ctScoreCovariate <- function(pieces, ptab, cov, G) {
  S <- pieces$S
  q <- ncol(G)
  k <- length(cov)
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  Stot <- colSums(S)
  U <- as.numeric(crossprod(G, S[, cov, drop = FALSE]))
  Idt <- matrix(0, q * k, ncol(S))
  Idd <- matrix(0, q * k, q * k)
  block <- function(a) (a - 1L) * q + seq_len(q)
  for (a in seq_len(k)) {
    Ha <- pieces$H[, match(cov[a], pieces$rows), , drop = TRUE]
    if (is.null(dim(Ha))) Ha <- matrix(Ha, nrow = nrow(S))
    Idt[block(a), ] <- -crossprod(G, Ha)
    for (b in seq_len(k)) {
      Idd[block(a), block(b)] <- -crossprod(G, G * Ha[, cov[b]])
    }
  }
  Ustar <- U - as.numeric(Idt %*% (Iinv %*% Stot))
  V <- Idd - Idt %*% Iinv %*% t(Idt)
  V <- (V + t(V)) / 2
  joint <- .ctScorePinv(V)
  epc <- as.numeric(joint$inv %*% Ustar)
  lm <- as.numeric(t(Ustar) %*% joint$inv %*% Ustar)
  rows <- lapply(seq_len(k), function(a) {
    b <- block(a)
    given <- .ctScorePinv(joint$inv[b, b, drop = FALSE])
    wald <- as.numeric(t(epc[b]) %*% given$inv %*% epc[b])
    alone <- .ctScorePinv(V[b, b, drop = FALSE])
    lmalone <- as.numeric(t(Ustar[b]) %*% alone$inv %*% Ustar[b])
    j <- cov[a]
    est <- if (q == 1L) epc[b] else NA_real_
    se <- if (q == 1L) sqrt(joint$inv[b, b]) else NA_real_
    data.frame(parameter = ptab$name[j], matrix = ptab$matrix[j],
      estimate = est, se = se, estimate_natural = ptab$slope[j] * est,
      df = given$rank,
      p = if (given$rank) stats::pchisq(wald, given$rank, lower.tail = FALSE) else NA_real_,
      p_alone = if (alone$rank) stats::pchisq(lmalone, alone$rank, lower.tail = FALSE) else NA_real_,
      stringsAsFactors = FALSE)
  })
  parameters <- do.call(rbind, rows)
  parameters$p_adj <- .ctScoreBH(parameters$p)
  list(joint = data.frame(covariate = attr(G, 'name'), form = attr(G, 'form'),
    basis = q, LM = lm, df = joint$rank,
    p = stats::pchisq(lm, joint$rank, lower.tail = FALSE),
    nmissing = attr(G, 'nmissing'), stringsAsFactors = FALSE),
    parameters = parameters[order(parameters$p), , drop = FALSE])
}

# The fitted model with a term `d * eta_k * eta_m` in each latent equation's
# CINT, for k <= m, every d a new PARS parameter with an identity transform.
# A CINT cell that is a free parameter moves into PARS, transform, variation
# and predictor effects and all, so its raw coordinate is unchanged; a fixed
# one keeps its value; an expression keeps its form.
.ctScoreStateModel <- function(model, order = 2L) {
  pars <- model$pars
  latent <- model$latentNames
  nl <- length(latent)
  # The monomials of degree 2 (and 3), each as its sorted state indices: the
  # first is the state the drift coefficient multiplies, the rest what it
  # varies with.
  monomials <- list()
  for (k in seq_len(nl)) for (m in k:nl) {
    monomials[[length(monomials) + 1L]] <- c(k, m)
    if (order >= 3L) for (l in m:nl) monomials[[length(monomials) + 1L]] <- c(k, m, l)
  }
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  ispars <- pars$matrix %in% 'PARS'
  nextpars <- max(c(0L, as.integer(pars$row[ispars])))
  template <- pars[which(pars$matrix %in% 'CINT')[1L], , drop = FALSE]
  parsrow <- function(name, from = NULL) {
    r <- if (is.null(from)) template else from
    r$matrix <- 'PARS'; r$row <- nextpars <<- nextpars + 1L; r$col <- 1L
    r$param <- name
    if (is.null(from)) {
      r$value <- NA; r$transform <- 'param'; r$indvarying <- FALSE
      if ('sdscale' %in% names(r)) r$sdscale <- 1
      for (e in effects) r[[e]] <- 'FALSE'
      for (l in levelflags) r[[l]] <- FALSE
    }
    r
  }
  terms <- list()
  added <- list()
  for (i in seq_len(nl)) {
    cell <- which(pars$matrix %in% 'CINT' & pars$row == i)
    if (length(cell) != 1L) stop("The model has no CINT cell for latent ", i,
      ".", call. = FALSE)
    names_i <- character()
    for (mono in monomials) {
      nm <- paste0('ctscoreterm', length(terms) + 1L)
      terms[[length(terms) + 1L]] <- data.frame(name = nm, equation = latent[i],
        i = i, k = mono[1L], by = paste(latent[mono[-1L]], collapse = '*'),
        term = paste(latent[mono], collapse = '*'),
        other = if (length(mono) == 2L && mono[1L] != mono[2L]) mono[2L] else NA_integer_,
        stringsAsFactors = FALSE)
      names_i <- c(names_i, paste(c(nm, latent[mono]), collapse = ' * '))
      added[[length(added) + 1L]] <- parsrow(nm)
    }
    sum_i <- paste(names_i, collapse = ' + ')
    current <- pars$param[cell]
    if (is.na(current) || !nzchar(current)) {
      value <- suppressWarnings(as.numeric(pars$value[cell]))
      expr <- if (is.na(value) || value == 0) sum_i else paste0(value, ' + ', sum_i)
    } else if (grepl('^[A-Za-z.][A-Za-z0-9._]*$', current) && !current %in% latent) {
      added[[length(added) + 1L]] <- parsrow(current, from = pars[cell, , drop = FALSE])
      expr <- paste0(current, ' + ', sum_i)
    } else {
      expr <- paste0('(', current, ') + ', sum_i)
    }
    pars$param[cell] <- expr
    pars$value[cell] <- NA
    pars$transform[cell] <- NA
    pars$indvarying[cell] <- FALSE
    for (e in effects) pars[[e]][cell] <- 'FALSE'
    for (l in levelflags) pars[[l]][cell] <- FALSE
  }
  model$pars <- rbind(pars, do.call(rbind, added))
  rownames(model$pars) <- NULL
  list(model = model, terms = do.call(rbind, terms))
}

# State dependence by the extended model: prepared as the fit was, evaluated
# with the fitted values and the terms at zero, which must reproduce the fit;
# then its gradient and its Hessian columns along the terms.
.ctScoreStateDependence <- function(fit, est, information = NULL, order = 2L,
  cores = 2L) {
  spec <- fit$model_spec
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) stop("This fit does not carry the model it was built from, ",
    "so the extended model cannot be prepared. Refit with this version of ctsem.",
    call. = FALSE)
  ext <- .ctScoreStateModel(base, order)
  args <- .ctFitReplayArgs(fit)
  args$model <- ext$model
  args$fit <- FALSE
  args$inits <- NULL
  args$datalong <- .ctFitLongData(fit)
  args$intoverpop <- .ctBackendIntOverPop(spec)
  args$verbose <- 0
  spec2 <- suppressMessages(suppressWarnings(do.call(ctFit, args)))
  nm1 <- as.character(.ctFitRawParNames(fit))
  nm2 <- as.character(.ctBackendRawParameterNames(list(model_spec = spec2),
    .ctBackendNpar(spec2)))
  lost <- setdiff(nm1, nm2)
  if (length(lost) || anyDuplicated(nm2)) stop("The extended model's parameters ",
    "do not line up with the fit's (", paste(lost, collapse = ', '), "). ",
    "State dependence cannot be tested for this model.", call. = FALSE)
  newnames <- setdiff(nm2, nm1)
  if (!setequal(newnames, ext$terms$name)) stop("The extended model added ",
    "parameters other than its terms: ", paste(setdiff(newnames, ext$terms$name),
      collapse = ', '), call. = FALSE)
  x2 <- setNames(numeric(length(nm2)), nm2)
  x2[nm1] <- est
  value1 <- ctJuliaEvaluate(fit, est, gradient = FALSE, cores = cores)$value
  ev <- ctJuliaEvaluate(spec2, x2, gradient = TRUE, cores = cores)
  if (!isTRUE(abs(ev$value - value1) <= 1e-6 * max(1, abs(value1)))) stop(
    "The extended model does not reproduce the fit at zero (log likelihood ",
    signif(ev$value, 8), " against ", signif(value1, 8), "), so its scores are ",
    "not the fit's. State dependence cannot be tested for this model.",
    call. = FALSE)
  grad <- setNames(as.numeric(ev$gradient), nm2)
  theta <- match(nm1, nm2)
  delta <- match(ext$terms$name, nm2)
  module <- .ctJuliaModule(spec$project)
  Hc <- .ctBackendWithMaxChunks(cores, as.matrix(.ctBackendJuliaValue(
    module$ctsem_hessian_columns(.ctJuliaObjective(spec2, hybrid = TRUE),
      .ctJuliaNumericVector(as.numeric(x2)), .ctJuliaVector(as.integer(delta))))))
  if (is.null(information)) information <- -.ctBackendHessian(fit, est)
  Iinv <- .ctScorePinv(information)$inv
  Idt <- -t(Hc[theta, , drop = FALSE])
  Idd <- -Hc[delta, , drop = FALSE]
  Idd <- (Idd + t(Idd)) / 2
  Ustar <- grad[delta] - as.numeric(Idt %*% (Iinv %*% grad[theta]))
  V <- Idd - Idt %*% Iinv %*% t(Idt)
  V <- (V + t(V)) / 2
  v <- diag(V)
  terms <- ext$terms
  latent <- .ctFitBaseModel(fit)$latentNames
  lm <- ifelse(v > 0, Ustar^2 / v, NA_real_)
  table <- data.frame(term = terms$term, equation = terms$equation,
    reads = paste0('DRIFT[', latent[terms$i], ',', latent[terms$k], '] varies with ',
      terms$by, ifelse(is.na(terms$other), '', paste0(', or DRIFT[',
        latent[terms$i], ',', latent[pmax(1L, terms$other)], '] with ',
        latent[terms$k]))),
    LM = lm, p = stats::pchisq(lm, 1, lower.tail = FALSE),
    estimate = ifelse(v > 0, Ustar / v, NA_real_),
    se = ifelse(v > 0, 1 / sqrt(v), NA_real_), stringsAsFactors = FALSE)
  table$p_adj <- .ctScoreBH(table$p)
  jointtest <- function(idx, label) {
    j <- .ctScorePinv(V[idx, idx, drop = FALSE])
    stat <- as.numeric(t(Ustar[idx]) %*% j$inv %*% Ustar[idx])
    data.frame(equation = label, LM = stat, df = j$rank,
      p = if (j$rank) stats::pchisq(stat, j$rank, lower.tail = FALSE) else NA_real_,
      stringsAsFactors = FALSE)
  }
  joint <- do.call(rbind, c(lapply(unique(terms$equation), function(e)
    jointtest(which(terms$equation == e), e)),
    list(jointtest(seq_along(delta), 'all'))))
  attr(table, 'joint') <- joint
  table <- table[order(table$p), , drop = FALSE]
  list(terms = table, joint = joint)
}

#' @export
print.ctScoreTest <- function(x, digits = 3, top = 10, ...) {
  fmt <- function(d) {
    d <- utils::head(d, top)
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], function(v) signif(v, digits))
    print(d, row.names = FALSE)
  }
  cat("Score tests at the fitted model (", x$info$route, " route, ",
    x$info$nsubjects, " subjects)\n", sep = "")
  if (!is.null(x$heterogeneity)) {
    cat("\nUnmodelled heterogeneity (one-sided; sd is the one-step estimate, raw scale):\n")
    fmt(x$heterogeneity[, c('parameter', 'z', 'p', 'p_adj', 'sd', 'sd_natural')])
  }
  for (c in x$covariates) {
    j <- c$joint
    cat("\nCovariate ", j$covariate, " (", j$form, "): joint LM = ",
      signif(j$LM, digits), ", df = ", j$df, ", p = ", signif(j$p, digits), "\n", sep = "")
    fmt(c$parameters[, c('parameter', 'estimate', 'se', 'estimate_natural', 'p',
      'p_adj', 'p_alone')])
  }
  if (!is.null(x$statedependence)) {
    cat("\nState dependence of the drift:\n")
    fmt(x$statedependence$joint)
    cat("\n")
    fmt(x$statedependence$terms[, c('term', 'equation', 'LM', 'p', 'p_adj',
      'estimate', 'se', 'reads')])
  }
  invisible(x)
}
