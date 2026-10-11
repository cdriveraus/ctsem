# Score tests of relaxations of a fitted model.
#
# Every relaxation tested here is an extra coordinate whose value zero is the
# fitted model: a random effect on a parameter that does not yet vary, a time
# independent covariate's effect on a parameter, a term making the drift depend
# on the latent state. The test needs only the score in that coordinate and the
# curvature blocks that involve it, evaluated at the fit -- no refit -- and
# projects the fitted parameters out of it (Neyman's C(alpha) form), so the
# test holds at an estimate that is near, not exactly at, the maximum (priors,
# a finish that stopped on gain). The design, the literature and the
# simulations behind the choices below are in
# CT-SEM/review/SCORETESTS-misfit-design-2026-10-09.md.
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
# The Hessian projects and sizes; it does not set the variance. A screen meets
# models that are wrong, and when one is wrong in any way the information
# equality fails, so a variance built from the Hessian mis-scales every other
# test: an unrelated covariate's joint test rejected 27 to 44 per cent of the
# time under a missing nonlinearity. The variance is the spread of each
# subject's projected contribution (Boos, 1992), and the p-value a studentised
# bootstrap of those contributions, since a variance score is a sum of skewed
# terms; the same unrelated test then rejected 2 per cent.
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
#' dynamics that depend on the latent state (a missing nonlinearity). Each test
#' is a score (Lagrange multiplier) test: the slope of the log likelihood in
#' the direction of the relaxation, at the fitted model, with the fitted
#' parameters projected out. Each also gives a one-step estimate of the
#' relaxation's size, for screening: refit with the relaxations that matter
#' before interpreting their sizes.
#'
#' @section Heterogeneity: For each parameter that is free, does not already
#'   vary over subjects, and is not an initial condition (T0MEANS and T0VAR,
#'   whose variation T0VAR already carries), the score test that its variance
#'   over subjects is zero. A variance cannot be negative, so the test is one
#'   sided: \code{z} large and positive says the subjects disagree about the
#'   parameter more than one shared value allows. The statistic is the
#'   information matrix test in that parameter, which is the score test for
#'   neglected heterogeneity (Chesher, 1984). It responds to any failure of
#'   the information equality in that coordinate, of which subject differences
#'   are the usual cause but not the only one: a subject-varying drift also
#'   moves the diffusion's test, because both change how much a subject
#'   varies, and a missing nonlinearity moves several. \code{sd} is the
#'   one-step estimate of the subjects' standard deviation on the raw scale,
#'   and \code{sd_natural} approximately on the parameter's own. To act on a
#'   result, set \code{indvarying = TRUE} for the parameter, refit, and test
#'   again. On a reduced-rank level (\code{poprank}), a parameter that varies
#'   only through the shared dimensions, without a variance of its own, is
#'   tested too: for variation outside those dimensions.
#'
#'   \code{dimension} tests every direction at once: whether subjects vary
#'   along some combination of the tested parameters -- a new dimension of
#'   person variation, or on a reduced-rank level one beyond its rank. Adding
#'   variance along a direction \code{u} changes a subject's log likelihood,
#'   at zero, by \code{(u's_i)^2/2 + u'H_iu/2}, so the per-parameter tests are
#'   the diagonal of one score matrix; its largest eigenvalue, scaled by each
#'   diagonal's spread, is the statistic, against a bootstrap over subjects,
#'   and its eigenvector is the direction, read as relative loadings (raw
#'   scale). It can detect variation spread over parameters that no single
#'   test singles out, and is less sensitive than the best single test when
#'   one parameter carries it.
#'
#' @section Covariates: For each covariate, the effects it would have on every
#'   tested parameter if added as a time independent predictor (on the raw
#'   scale, as a predictor effect acts), including on the population standard
#'   deviations of the parameters that vary: a covariate moderating
#'   heterogeneity. The \code{joint} test of all of a covariate's effects is the
#'   screen. Per parameter, \code{estimate} is the one-step effect per unit of
#'   the covariate with the fitted parameters held at their estimates, jointly
#'   over the parameters, so a covariate that shifts one parameter is not read
#'   off another correlated with it; \code{estimate_natural} puts it
#'   approximately on the parameter's own scale, and \code{p} tests it. Holding
#'   the fit's parameters is what keeps the size near a refit's: a covariate
#'   the model lacks is absorbed by the parameters it has -- a random effect's
#'   standard deviation grows to take in a covariate that shifts its mean --
#'   and letting them move in a one-step extrapolates the likelihood across
#'   that whole change. A covariate enters linearly (centred), as treatment
#'   indicators for a factor, or as a natural spline with \code{df} degrees of
#'   freedom, for which only \code{p} is given; subjects missing it contribute
#'   at the mean. To act on a result, add the covariate to the data and to
#'   \code{TIpredNames} and refit.
#'
#'   \code{form = 'threshold'} tests a step at an unknown cut instead: the
#'   joint test at each cut from the 15th to the 85th percentile, the largest
#'   of them against a bootstrap over subjects that takes the search over cuts
#'   into account (Andrews, 1993) -- the split a score-guided tree makes. The
#'   effects reported are those of being above the best cut, \code{cut} in the
#'   joint row.
#'
#'   On a reduced-rank level the spread of a varying parameter has no
#'   \code{popsd_} coordinate, so a covariate's effect on it is tested along
#'   the direction that scales that parameter's row of the loadings and its own
#'   sd together, reported as \code{scale(<parameter>)}: the effect on its log
#'   standard deviation per unit of the covariate, the correlations held.
#'
#' @section State dependence: For each latent equation and each pair of
#'   latent states \code{k <= m}, the score test of a term \code{d * eta_k *
#'   eta_m} in that equation's drift: a drift coefficient that varies linearly
#'   with a state, the second-order expansion of a drift that depends on the
#'   state (Luukkonen, Saikkonen & Terasvirta, 1988). A term \code{eta_k *
#'   eta_m} in the equation of \code{eta_i} is \code{DRIFT[i,k]} varying with
#'   \code{eta_m} -- and, when \code{k != m}, equally \code{DRIFT[i,m]} varying
#'   with \code{eta_k}; at this order the data cannot tell those apart, and
#'   \code{reads} says both. \code{estimate} is the one-step change in that
#'   drift coefficient per unit of the state. It is \code{NA} when the fitted model's
#'   curvature cannot size the term, which happens when the departure is
#'   strong; the test itself is unaffected. \code{order = 3} adds the terms
#'   of degree three. The terms are evaluated by the model's own filter at the
#'   state predicted from earlier observations, which keeps the test valid; a
#'   state estimate that has seen the observation being scored -- the updated
#'   or smoothed state -- would not. To act on a result, write the dependence
#'   into the drift cell (\code{DRIFT[1,1] = 'dr11 + d * eta2'} with
#'   \code{PARS = c('dr11|-log1p_exp(-param)', 'd')}) and refit.
#'
#' @section Fixed cells: For each cell the model fixes at zero that could
#'   sensibly be freed -- a \code{DRIFT} cross effect, a \code{DIFFUSION} or
#'   \code{MANIFESTVAR} covariance, a \code{LAMBDA} cross loading, a
#'   \code{TDPREDEFFECT} -- the score test of freeing it: the modification
#'   index (Sorbom, 1989), with \code{estimate} the expected parameter change.
#'   \code{CINT} and \code{MANIFESTMEANS} are not offered, since with the
#'   other free either is identified only through the initial transient.
#'   State dependence and fixed cells are tested on one extended model, so
#'   asking for both costs one preparation.
#'
#' @section Stability: \code{tests = 'instability'}, not run unless asked
#'   for, tests whether each parameter is stable across the data: ordered by
#'   time, and by each latent's state as predicted from earlier observations.
#'   It uses each data row's score -- the derivative of its one-step log
#'   likelihood -- whose partial sums, under the model, behave as a Brownian
#'   bridge in any order fixed in advance (Nyblom, 1989). A parameter that
#'   changes over the study, or differs at high and low levels of a state,
#'   bends them. Ordered by a state, that is a test of state dependence of
#'   unknown form, beside the polynomial terms of \code{statedependence}. The
#'   estimation effect is removed with the information accumulated along the
#'   ordering, and the p-values are from a wild bootstrap over rows. It costs a
#'   dual pass of the filter per eight parameters, and on the Laplace route one
#'   per truncated record of every subject.
#'
#' @section Alike relaxations: Two relaxations whose per-subject projected
#'   contributions are nearly collinear (correlation of 0.9 or more) would
#'   each explain the same signal, and the data cannot tell them apart at this
#'   order (MacCallum et al., 1993). They share a letter in \code{class};
#'   \code{NA} means alike to none. Refit with the most plausible, then test
#'   again.
#'
#' @section Variance and p-values: The fitted model's curvature projects the
#'   fitted parameters out of each test and sets the one-step sizes; the
#'   variance of each test is, by default (\code{variance = 'robust'}), the
#'   spread over subjects of their projected contributions, and its p-value a
#'   bootstrap over subjects of the studentised statistic (\code{nboot}
#'   resamples, from a fixed stream, so a call is reproducible and leaves the
#'   session's random numbers alone). That keeps each test valid when the model
#'   is wrong in some other way, which is when a screen is used, and corrects
#'   for the skewness of a variance score. \code{variance = 'model'} uses the
#'   curvature for the variance and chi-square and normal references instead;
#'   it is what is used below 20 subjects, where a spread over subjects is too
#'   thin to estimate.
#'
#' @section Reading the results: The tests are screens. A badly misspecified
#'   model lights up many at once, and their one-step sizes are local
#'   approximations -- especially for drift parameters, which act through the
#'   matrix exponential. Relax what stands out most, refit, and test again.
#'   Within each family \code{p_adj} is the Benjamini-Hochberg adjustment over
#'   the tests shown. The tests use the fit's own objective, so they apply to
#'   the route it was fitted on (\code{intoverpop}); they are made at an
#'   optimum, and refuse a sampled fit. Heterogeneity and covariate tests need
#'   subjects as the independent units, so a model with a level above the
#'   subject can only be tested for state dependence. A family that cannot be
#'   tested on a model says why, as a warning and in \code{info$errors}, and
#'   the other families' results stand. With \code{priors = TRUE} the new terms
#'   of the extended model would carry the normal(0, 1) prior every parameter
#'   has; they are tested on the likelihood, without it.
#'
#' @param fit A fit from \code{\link{ctFit}} with \code{backend = 'julia'},
#'   optimised.
#' @param covariates Optional data frame of subject-level covariates to test:
#'   a column named as the fit's subject id (or \code{id}) and one column per
#'   covariate. Long data with repeated rows per subject is accepted when each
#'   covariate is constant within subject. Without an id column, one row per
#'   subject in the fit's subject order.
#' @param tests Which families to run: any of \code{'heterogeneity'},
#'   \code{'covariates'}, \code{'statedependence'}, \code{'fixed'},
#'   \code{'instability'}. By default all but \code{'instability'}; covariates
#'   are tested only when \code{covariates} is given.
#' @param form How a covariate enters: \code{'auto'} (indicators for a factor,
#'   character or logical covariate; linear otherwise), \code{'linear'},
#'   \code{'spline'}, \code{'factor'} or \code{'threshold'} (a step at an
#'   unknown cut). One value, or one per covariate.
#' @param df Degrees of freedom of the spline basis.
#' @param parameters Optional parameter names to restrict the heterogeneity and
#'   covariate tests to.
#' @param order Highest degree of the state-dependence terms: 2, each drift
#'   coefficient varying linearly with each state; 3 adds variation with
#'   products and squares of states, which is what a symmetric dependence --
#'   regulation that strengthens at both extremes, say -- shows up as, and
#'   which degree 2 cannot see.
#' @param variance \code{'robust'} (the default) or \code{'model'}; see the
#'   section on variance.
#' @param nboot Bootstrap resamples for the robust p-values.
#' @param cores Engine cores for the evaluations.
#' @param verbose Report progress.
#' @return An object of class \code{ctScoreTest}: a list with
#'   \code{heterogeneity} (a data frame, one row per parameter),
#'   \code{dimension} (the largest eigenvalue, its p-value and the direction,
#'   when two or more parameters are tested),
#'   \code{covariates} (one element per covariate, each a list of its
#'   \code{joint} test and a data frame \code{parameters}),
#'   \code{statedependence} (a list of \code{terms}, a data frame, and
#'   \code{joint}, the tests per equation and overall), \code{fixed} (a list
#'   of \code{cells} and their \code{joint} test), \code{instability} (a list
#'   of \code{parameters}, one row per ordering and parameter, and
#'   \code{joint}, per ordering), and \code{info} (the
#'   route, the number of subjects and parameters, the variance used, and
#'   \code{errors}, why any family was not tested). Tables carry \code{class}
#'   where relaxations are alike.
#' @references Boos, D. D. (1992). On generalized score tests. The American
#'   Statistician, 46, 327-333.
#'
#'   Chesher, A. (1984). Testing for neglected heterogeneity. Econometrica, 52,
#'   865-872.
#'
#'   Andrews, D. W. K. (1993). Tests for parameter instability and structural
#'   change with unknown change point. Econometrica, 61, 821-856.
#'
#'   Nyblom, J. (1989). Testing for the constancy of parameters over time.
#'   Journal of the American Statistical Association, 84, 223-230.
#'
#'   Hall, D. B., & Praestgaard, J. T. (2001). Order-restricted score tests for
#'   homogeneity in generalised linear and nonlinear mixed models.
#'   Biometrika, 88, 739-751.
#'
#'   MacCallum, R. C., Wegener, D. T., Uchino, B. N., & Fabrigar, L. R. (1993).
#'   The problem of equivalent models in applications of covariance structure
#'   analysis. Psychological Bulletin, 114, 185-199.
#'
#'   Luukkonen, R., Saikkonen, P., & Terasvirta, T. (1988). Testing linearity
#'   against smooth transition autoregressive models. Biometrika, 75, 491-499.
#'
#'   Sorbom, D. (1989). Model modification. Psychometrika, 54, 371-384.
#' @seealso \code{vignette('score-tests')} for three worked examples.
#' @examples
#' \dontrun{
#' fit <- ctFit(dat, model, backend = 'julia')
#' st <- ctScoreTest(fit, covariates = subjectdata)
#' st
#' st$heterogeneity
#' }
#' @export
ctScoreTest <- function(fit, covariates = NULL,
  tests = c('heterogeneity', 'covariates', 'statedependence', 'fixed', 'instability'),
  form = c('auto', 'linear', 'spline', 'factor', 'threshold'), df = 3, parameters = NULL,
  order = 2L, variance = c('robust', 'model'), nboot = 1999L, cores = 2L,
  verbose = TRUE) {
  if (!.ctFitIsJulia(fit)) stop("ctScoreTest() needs a fit from ",
    "ctFit(backend = 'julia'). Refit on that backend to use it.", call. = FALSE)
  if (!is.null(fit$sample)) stop("Score tests are made at an optimum, and this ",
    "fit was sampled. Fit with optimize = TRUE to test it.", call. = FALSE)
  if (isFALSE(fit$args$resolved$intoverstates)) stop("intoverstates = FALSE ",
    "targets the joint density of parameters and states, whose maximum is not ",
    "an estimate of anything (see ?ctFit); there is no fitted model to test.",
    call. = FALSE)
  # Before anything reaches Julia: a session this call starts is as wide as
  # `cores`, as ctFit's is. Without it a saved fit tested in a fresh R session
  # (a cluster job) ran every pass at getOption('mc.cores', 2) threads.
  .ctBackendResolveThreads(cores)
  testsgiven <- !missing(tests)
  tests <- match.arg(tests, several.ok = TRUE)
  # Stability needs a dual pass of the filter per eight parameters (on the
  # Laplace route one per truncated record), so it is run when asked for.
  if (!testsgiven) tests <- setdiff(tests, 'instability')
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
  variance <- match.arg(variance)
  if (variance == 'robust' && .ctFitNsubjects(fit) < 20L) {
    say("Fewer than 20 subjects: using the model-based variance.")
    variance <- 'model'
  }
  # The bootstrap draws from its own stream, so a call gives the same answer
  # every time and leaves the caller's random numbers where they were.
  if (exists('.Random.seed', envir = globalenv(), inherits = FALSE)) {
    oldseed <- get('.Random.seed', envir = globalenv(), inherits = FALSE)
    on.exit(assign('.Random.seed', oldseed, envir = globalenv()), add = TRUE)
  } else on.exit(rm('.Random.seed', envir = globalenv()), add = TRUE)
  set.seed(20261009L)

  spec <- fit$model_spec
  est <- .ctFitRawEstimate(fit)
  ptab <- .ctScoreParameterTable(fit)
  out <- list(heterogeneity = NULL, covariates = NULL, statedependence = NULL,
    fixed = NULL, info = list(route = .ctBackendIntOverPop(spec),
      nsubjects = .ctFitNsubjects(fit), npar = length(est)))
  out$info$variance <- variance
  out$info$nboot <- nboot
  # A family that cannot be tested on this model says why, in the result and
  # as a warning, and the others stand: the pilot's state-dependence refusal
  # once discarded a heterogeneity family that had taken twenty minutes.
  attempt <- function(label, expr) {
    ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
      out$info$errors[[label]] <<- conditionMessage(e)
      warning(label, ": ", conditionMessage(e), call. = FALSE)
      FALSE
    })
    invisible(ok)
  }

  persubject <- intersect(tests, c('heterogeneity', 'covariates'))
  pieces <- NULL
  # A reduced-rank level (poprank) carries its variances as loadings and own
  # sds, so a parameter already varying can still lack a dimension of
  # variation, and a covariate's effect on its spread has no popsd coordinate
  # to act on: it acts on the scale of the parameter's row of the loadings.
  reduced <- any(grepl('^poploading_', ptab$name))
  if (length(persubject)) {
    het <- .ctScoreSelect(ptab, parameters, 'heterogeneity', reduced)
    cov <- .ctScoreSelect(ptab, parameters, 'covariates', reduced)
    if (!'heterogeneity' %in% tests) het <- integer()
    if (!'covariates' %in% tests) cov <- integer()
    # The covariates and their bases are checked before any pass is paid for:
    # a wrong id column or a non-numeric threshold is the caller's to fix, and
    # stops the call, where a family the model cannot support does not.
    bases <- NULL
    if (length(cov)) {
      values <- .ctScoreCovariateValues(fit, covariates)
      forms <- rep_len(form, length(values))
      bases <- Map(function(z, f, name) {
        if (identical(f, 'threshold')) {
          if (is.factor(z) || is.character(z) || is.logical(z)) stop("Covariate '", name,
            "' is not numeric; a threshold needs an ordered covariate.", call. = FALSE)
          list(threshold = z)
        } else .ctScoreBasis(z, f, df, name)
      }, values, forms, names(values))
    }
    dirs <- if (reduced && length(cov)) .ctScoreScaleDirections(est, ptab, parameters) else NULL
    say("Score contributions: ", 2L * (length(union(het, cov)) + NCOL(dirs)) + 1L,
      " per-subject gradient passes")
    pieces <- .ctScorePieces(fit, est, rows = union(het, cov), second = het,
      dirs = dirs, cores = cores)
    if (nrow(pieces$S) != .ctFitNsubjects(fit)) stop("This model has a level ",
      "above the subject, so its independent units are not subjects. ",
      "Heterogeneity and covariate tests need subjects; use ",
      "tests = 'statedependence'.", call. = FALSE)
    if (length(het)) attempt('heterogeneity', {
      out$heterogeneity <- .ctScoreHeterogeneity(pieces, ptab, het, variance, nboot)
      if (length(het) > 1L) out$dimension <- .ctScoreDimension(pieces, ptab, het, nboot)
    })
    if (length(cov)) attempt('covariates', {
      items <- .ctScoreItems(pieces, ptab, cov)
      out$covariates <- Map(function(G, name) if (is.list(G) && !is.null(G$threshold))
        .ctScoreThreshold(pieces, items, G$threshold, name, variance, nboot) else
        .ctScoreCovariate(pieces, items, G, variance, nboot), bases, names(bases))
    })
  }
  extended <- intersect(tests, c('statedependence', 'fixed'))
  if (length(extended)) attempt(paste(extended, collapse = ' and '), {
    say("Preparing the model extended at zero (", paste(extended, collapse = ', '), ")")
    information <- if (!is.null(pieces)) -pieces$Hsum else NULL
    ext <- .ctScoreExtended(fit, est, extended, information, order = order,
      variance = variance, nboot = nboot, cores = cores)
    out$statedependence <- ext$statedependence
    out$fixed <- ext$fixed
  })
  if ('instability' %in% tests) attempt('instability', {
    say("Row scores for the stability tests")
    out$instability <- .ctScoreInstability(fit, est, ptab, parameters, nboot, cores)
  })
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

# The raw positions a family tests: model parameters for both, minus the
# initial conditions and those already varying for heterogeneity (on a
# reduced-rank level, those with a variance of their own); plus the population
# scales for covariates. `parameters` narrows either.
.ctScoreSelect <- function(ptab, parameters, family, reduced = FALSE) {
  # On a reduced-rank level a varying parameter with no variance of its own
  # (no popunique_) has none outside the shared dimensions: one of its own is
  # a boundary, testable as for a parameter that does not vary. A parameter
  # with its own variance already has that direction free, so a test there
  # asks nothing the fit did not estimate.
  own <- paste0('popunique_', ptab$name) %in% ptab$name
  keep <- if (family == 'heterogeneity') {
    ptab$kind == 'model' & (!ptab$indvarying | (reduced & !own)) &
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
# gradients along the tested parameters only (2 * length(rows) + 1 passes): a
# subject's Hessian is symmetric, so its column along a tested parameter is the
# row the tests read. The same passes give, for the heterogeneity parameters,
# the second derivative of every score component along that parameter, which
# is the third derivative the heterogeneity test's projection needs:
# d H_i[j,j] / d theta_k = d^2 s_ik / d theta_j^2. A step of 1e-3 on the raw
# scale balances truncation (h^2) against the rounding of the second
# difference (eps / h^2). The full information is the fit's own Hessian, which
# the optimiser's finish formed already; differencing along every parameter
# for it cost 2p + 1 passes, 705 on a 352-parameter fit.
.ctScorePieces <- function(fit, est, rows, second, dirs = NULL, step = 1e-3,
  cores = 2L) {
  p <- length(est)
  if (!length(rows)) stop("No parameters to test.", call. = FALSE)
  Hdir <- NULL
  .ctBackendWithMaxChunks(cores, {
    S <- .ctBackendScoreMatrix(fit, est)
    n <- nrow(S)
    H <- array(0, c(n, length(rows), p))
    D2 <- array(0, c(n, p, length(second)))
    for (a in seq_along(rows)) {
      k <- rows[a]
      h <- step * max(1, abs(est[k]))
      e <- est; e[k] <- est[k] + h
      plus <- .ctBackendScoreMatrix(fit, e)
      e[k] <- est[k] - h
      minus <- .ctBackendScoreMatrix(fit, e)
      H[, a, ] <- (plus - minus) / (2 * h)
      j <- match(k, second)
      if (!is.na(j)) D2[, , j] <- (plus - 2 * S + minus) / h^2
    }
    # Along each direction, the same difference gives H_i w.
    if (!is.null(dirs)) {
      Hdir <- array(0, c(n, ncol(dirs), p))
      for (d in seq_len(ncol(dirs))) {
        w <- dirs[, d]
        h <- step * max(1, max(abs(est[w != 0]))) / max(abs(w))
        plus <- .ctBackendScoreMatrix(fit, est + h * w)
        minus <- .ctBackendScoreMatrix(fit, est - h * w)
        Hdir[, d, ] <- (plus - minus) / (2 * h)
      }
    }
    Hsum <- .ctBackendHessian(fit, est)
  })
  Hsum <- (Hsum + t(Hsum)) / 2
  list(S = S, H = H, rows = rows, D2 = D2, second = second, Hsum = Hsum,
    dirs = dirs, Hdir = Hdir)
}

# The raw-coordinate direction that scales one varying parameter's spread on a
# reduced-rank level: its row of the loadings (`poploading_<par>_dim<q>`, raw
# times the parameter's sd scale) and its own sd if it has one
# (`popunique_<par>`, log1p_exp(2 raw - 1) times the scale; see
# `_laplace_poploading`). Multiplying that row by exp(g) moves each loading's
# raw value by g times itself and the own sd's by g f / f'. A covariate's
# effect on g is its effect on the parameter's log sd, the correlations held.
.ctScoreScaleDirections <- function(est, ptab, parameters = NULL) {
  nm <- ptab$name
  load <- grepl('^poploading_.*_dim[0-9]+$', nm)
  own <- grepl('^popunique_', nm) & !grepl('\\.', nm)
  pars <- unique(c(sub('_dim[0-9]+$', '', sub('^poploading_', '', nm[load])),
    sub('^popunique_', '', nm[own])))
  if (!is.null(parameters)) pars <- intersect(pars, parameters)
  W <- matrix(0, length(est), 0)
  for (par in pars) {
    w <- numeric(length(est))
    l <- which(load & startsWith(nm, paste0('poploading_', par, '_dim')) &
      grepl('^[0-9]+$', substring(nm, nchar(paste0('poploading_', par, '_dim')) + 1L)))
    w[l] <- est[l]
    u <- which(own & nm == paste0('popunique_', par))
    if (length(u)) {
      x <- 2 * est[u] - 1
      w[u] <- log1p(exp(x)) / (2 * stats::plogis(x))
    }
    if (any(w != 0)) W <- cbind(W, matrix(w, ncol = 1, dimnames = list(NULL, paste0('scale(', par, ')'))))
  }
  if (!ncol(W)) NULL else W
}

# What a covariate's effects act on, each with its score, its Hessian row and
# its direction in raw coordinates: the tested parameters, and on a
# reduced-rank level the scale of each varying parameter's spread.
.ctScoreItems <- function(pieces, ptab, cov) {
  p <- ncol(pieces$S)
  items <- lapply(cov, function(j) list(name = ptab$name[j], matrix = ptab$matrix[j],
    slope = ptab$slope[j], s = pieces$S[, j], H = .ctScoreHrow(pieces, j),
    dir = replace(numeric(p), j, 1)))
  if (!is.null(pieces$dirs)) for (d in seq_len(ncol(pieces$dirs))) {
    H <- pieces$Hdir[, d, , drop = TRUE]
    if (is.null(dim(H))) H <- matrix(H, nrow = nrow(pieces$S))
    items[[length(items) + 1L]] <- list(name = colnames(pieces$dirs)[d],
      matrix = 'POPCOV', slope = 1, s = as.numeric(pieces$S %*% pieces$dirs[, d]),
      H = H, dir = pieces$dirs[, d])
  }
  items
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
# subjects, projected on the fitted parameters (C(alpha)): each subject's
# contribution plus D' I^-1 s_i, with D the derivative of the statistic in
# theta. Its variance, `robust`, is the spread of those contributions, which
# stays right when the model is wrong in some other way -- the case a screen
# meets most -- where the information equality the model-based variance
# leans on fails. `model` is sum_i H_i,jj^2 / 2 less D' I^-1 D: exact for a
# Gaussian score, as a parameter of a linear Gaussian model's mean has, and
# the information for the variance in either case, which sets the one-step
# estimate of the standard deviation.
.ctScoreHeterogeneity <- function(pieces, ptab, het, variance = 'robust',
  nboot = 1999L) {
  S <- pieces$S
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  Stot <- colSums(S)
  rows <- lapply(seq_along(het), function(a) {
    j <- het[a]
    s <- S[, j]
    Hj <- .ctScoreHrow(pieces, j)
    hjj <- Hj[, j]
    # D = dT/dtheta: d(s_ij^2/2)/dtheta = s_ij H_i[j, ] and
    # d(H_i,jj/2)/dtheta_k = d^2 s_ik / dtheta_j^2 / 2.
    D <- colSums(s * Hj) + 0.5 * colSums(pieces$D2[, , a, drop = TRUE])
    w <- as.numeric(Iinv %*% D)
    # T(theta) + D' I^-1 S(theta) is T at the true theta to first order,
    # wherever it is evaluated. Plus: D is a derivative and I the negative
    # Hessian; the covariate test's `- Idt I^-1 S` is the same with
    # Idt = -dU/dtheta.
    contribution <- 0.5 * (s^2 + hjj) + as.numeric(S %*% w)
    stat <- sum(contribution)
    information <- sum(0.5 * hjj^2) - sum(D * w)
    V <- if (variance == 'robust') .ctScoreSpread(matrix(contribution))[1, 1] else
      information
    z <- if (V > 0) stat / sqrt(V) else NA_real_
    p <- if (variance == 'robust') .ctScoreBoot(matrix(contribution), nboot,
      onesided = TRUE) else stats::pnorm(z, lower.tail = FALSE)
    # The one-step variance divides by the information for it. Less what
    # estimating the other parameters takes, when that leaves it positive;
    # the third derivatives that projection needs are differenced twice and
    # can overshoot, and the plain sum is the size to first order anyway.
    if (information <= 0) information <- sum(0.5 * hjj^2)
    sdraw <- if (information > 0) sqrt(max(0, stat / information)) else NA_real_
    row <- data.frame(parameter = ptab$name[j], matrix = ptab$matrix[j],
      row = ptab$row[j], col = ptab$col[j], z = z, p = p, sd = sdraw,
      sd_natural = abs(ptab$slope[j]) * sdraw, stringsAsFactors = FALSE)
    attr(row, 'contribution') <- contribution
    row
  })
  contributions <- do.call(cbind, lapply(rows, attr, 'contribution'))
  out <- do.call(rbind, rows)
  out$p_adj <- .ctScoreBH(out$p)
  out$class <- .ctScoreClasses(contributions)
  out[order(out$p, -out$z), , drop = FALSE]
}

# Whether the subjects vary along a direction no single parameter's test
# isolates -- a new dimension of person variation over several parameters, or
# on a reduced-rank level one more than its rank. Adding variance tau^2 along
# a raw direction u changes subject i's log likelihood, at tau^2 = 0, by
# (u' s_i)^2 / 2 + u' H_i u / 2: the expectation over the subject's posterior
# of the conditional score's square and Hessian is the marginal's square plus
# its Hessian. So T = sum_i (s_i s_i' + H_i) / 2 over the tested parameters is
# the score matrix of every such direction; its diagonal is the
# per-parameter test. The off-diagonal projections omit the third
# derivatives the diagonal has, which can only widen the spread and so keeps
# the test conservative. The statistic is the largest eigenvalue of T scaled
# by each diagonal's spread (invariant to each parameter's units), against a
# multiplier bootstrap over subjects; its eigenvector is the direction, read
# as loadings (Hall & Praestgaard, 2001; Zhu & Zhang, 2006).
.ctScoreDimension <- function(pieces, ptab, het, nboot = 1999L) {
  S <- pieces$S
  n <- nrow(S); K <- length(het)
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  Hr <- lapply(het, function(j) .ctScoreHrow(pieces, j))
  C <- array(0, c(n, K, K))
  for (a in seq_len(K)) for (b in a:K) {
    ja <- het[a]; jb <- het[b]
    base <- 0.5 * (S[, ja] * S[, jb] + Hr[[a]][, jb])
    D <- 0.5 * (colSums(S[, ja] * Hr[[b]]) + colSums(S[, jb] * Hr[[a]]))
    if (a == b) D <- colSums(S[, ja] * Hr[[a]]) +
      0.5 * colSums(pieces$D2[, , match(ja, pieces$second), drop = TRUE])
    C[, a, b] <- C[, b, a] <- base + as.numeric(S %*% (Iinv %*% D))
  }
  T <- apply(C, c(2, 3), sum)
  dsd <- vapply(seq_len(K), function(a) sqrt(.ctScoreSpread(matrix(C[, a, a]))[1, 1]), 0)
  if (any(!is.finite(dsd) | dsd <= 0)) return(NULL)
  e <- eigen(T / sqrt(outer(dsd, dsd)), symmetric = TRUE)
  lambda <- e$values[1]
  # The same statistic from any rows of the flattened contributions: each
  # diagonal's spread recomputed, so the reference is studentised as lambda is.
  diagcols <- (seq_len(K) - 1L) * K + seq_len(K)
  largest <- function(X) {
    d <- sqrt(pmax(diag(.ctScoreSpread(X[, diagcols, drop = FALSE])), 0))
    if (any(!is.finite(d) | d <= 0)) return(NA_real_)
    eigen(matrix(colSums(X), K) / sqrt(outer(d, d)), symmetric = TRUE,
      only.values = TRUE)$values[1]
  }
  u <- e$vectors[, 1] / sqrt(dsd)
  u <- u / u[which.max(abs(u))]
  list(lambda = lambda, p = .ctScoreBoot(matrix(C, n, K * K), nboot, statistic = largest),
    direction = stats::setNames(u, ptab$name[het])[order(-abs(u))])
}

# Subject i's Hessian row for raw parameter j, n x p.
.ctScoreHrow <- function(pieces, j) {
  Hj <- pieces$H[, match(j, pieces$rows), , drop = TRUE]
  if (is.null(dim(Hj))) Hj <- matrix(Hj, nrow = nrow(pieces$S))
  Hj
}

# The variance of a sum of independent units' contributions (rows of C),
# estimated from their spread: n / (n - 1) times the centred cross product.
.ctScoreSpread <- function(C) {
  C <- as.matrix(C)
  n <- nrow(C)
  centred <- sweep(C, 2, colMeans(C))
  crossprod(centred) * n / (n - 1)
}

# The p-value of U' V^-1 U, U = colSums(C) and V its robust variance, by a
# studentised bootstrap over the units: rows of C centred (the null) and
# resampled, the statistic recomputed with its own variance each time. The
# chi-square (or normal) reference assumes U is Gaussian, and a variance
# score is a sum of chi-square-like terms: with 100 subjects the robust
# one-sided z rejected a true model at 2 per cent where 5 was asked, and the
# joint covariate test at 7.5 to 10.5. Studentising and resampling corrects
# the skewness to second order, and costs milliseconds: nothing is refitted.
# `onesided` takes a single column and rejects for large U. `statistic`
# replaces U' V^-1 U with another function of the contributions, recomputed the
# same way. `C` may be a list of contribution matrices over the same units, for
# the largest of their statistics -- a search, as over a threshold's cuts --
# against the largest in each resample.
#
# The draws must be studentised as the observed statistic is. A Gaussian
# multiplier with the observed variance held fixed is not, and errs either
# way: for the threshold search, whose sums are roughly symmetric, the fixed
# variance made the reference too narrow (21 per cent at .05 under the
# null); for the dimension test, whose diagonal sums are right-skewed and so
# left-skewed once studentised, too wide (0 of 100, and 2 per cent power for a
# missing dimension).
.ctScoreBoot <- function(C, nboot = 1999L, onesided = FALSE, statistic = NULL) {
  Cs <- if (is.list(C) && !is.data.frame(C)) lapply(C, as.matrix) else list(as.matrix(C))
  n <- nrow(Cs[[1]])
  stat <- if (!is.null(statistic)) statistic else function(X) {
    u <- colSums(X)
    v <- .ctScoreSpread(X)
    if (onesided) return(u / sqrt(v[1, 1]))
    pv <- .ctScorePinv(v)
    as.numeric(t(u) %*% pv$inv %*% u)
  }
  sup <- function(Xs) max(vapply(Xs, stat, numeric(1)))
  observed <- sup(Cs)
  if (!is.finite(observed)) return(NA_real_)
  nulls <- lapply(Cs, function(X) sweep(X, 2, colMeans(X)))
  draws <- vapply(seq_len(nboot), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    sup(lapply(nulls, function(X) X[idx, , drop = FALSE]))
  }, numeric(1))
  (1 + sum(draws >= observed, na.rm = TRUE)) / (1 + sum(is.finite(draws)))
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
# basis G: each subject's contribution g_i s_ij, projected on the fitted
# parameters with the per-subject Hessians weighted by G. The variance is the
# contributions' spread (`robust`) or the projected information (`model`).
#
# The effects are estimated with the fitted parameters held (Idd^-1 U*), not
# free (V^-1 U*, the classical expected parameter change). A covariate the
# model lacks is absorbed by the parameters it has -- a random effect's sd
# grows to take in a covariate that shifts its mean -- and the free one-step
# extrapolates the likelihood's curvature across that whole change. On
# simulated data where a covariate shifted a random manifest mean by 0.8:
# held 0.81, refit 0.81, free 3.0, and Newton from the fit diverged; over 100
# such datasets the free one also attributed effects to parameters the
# covariate did not touch (T0MEANS 50 per cent, a population sd 58). Holding
# them is bigIRT's post-hoc adjustment (its items held, the ability side
# estimated), here joint over all of the covariate's effects.
.ctScoreCovariate <- function(pieces, items, G, variance = 'robust',
  nboot = 1999L) {
  S <- pieces$S
  q <- ncol(G)
  k <- length(items)
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  block <- function(a) (a - 1L) * q + seq_len(q)
  # Contributions g_ia s_ij, ordered by item then basis column.
  C <- do.call(cbind, lapply(items, function(it) G * it$s))
  U <- colSums(C)
  Idt <- matrix(0, q * k, ncol(S))
  Idd <- matrix(0, q * k, q * k)
  for (a in seq_len(k)) {
    Ha <- items[[a]]$H
    Idt[block(a), ] <- -crossprod(G, Ha)
    for (b in seq_len(k)) Idd[block(a), block(b)] <-
      -crossprod(G, G * as.numeric(Ha %*% items[[b]]$dir))
  }
  Idd <- (Idd + t(Idd)) / 2
  A <- Idt %*% Iinv
  Cstar <- C - S %*% t(A)
  Ustar <- colSums(Cstar)
  Vmodel <- Idd - A %*% t(Idt)
  Vmodel <- (Vmodel + t(Vmodel)) / 2
  V <- if (variance == 'robust') .ctScoreSpread(Cstar) else Vmodel
  joint <- .ctScorePinv(V)
  lm <- as.numeric(t(Ustar) %*% joint$inv %*% Ustar)
  pjoint <- if (variance == 'robust') .ctScoreBoot(Cstar, nboot) else
    stats::pchisq(lm, joint$rank, lower.tail = FALSE)
  heldinv <- .ctScorePinv(Idd)$inv
  held <- as.numeric(heldinv %*% Ustar)
  heldcov <- heldinv %*% V %*% heldinv
  # Each subject's contribution to the held estimates, for their bootstrap:
  # summed they are the estimates, and their spread is `heldcov`.
  heldC <- Cstar %*% t(heldinv)
  wald <- function(b, est, covariance) {
    if (variance == 'robust') return(.ctScoreBoot(heldC[, b, drop = FALSE], nboot))
    w <- .ctScorePinv(covariance[b, b, drop = FALSE])
    if (!w$rank) return(NA_real_)
    stats::pchisq(as.numeric(t(est[b]) %*% w$inv %*% est[b]), w$rank, lower.tail = FALSE)
  }
  rows <- lapply(seq_len(k), function(a) {
    b <- block(a)
    it <- items[[a]]
    one <- function(x) if (q == 1L) x else NA_real_
    data.frame(parameter = it$name, matrix = it$matrix,
      estimate = one(held[b]), se = one(sqrt(heldcov[b, b])),
      estimate_natural = it$slope * one(held[b]),
      p = wald(b, held, heldcov), stringsAsFactors = FALSE)
  })
  parameters <- do.call(rbind, rows)
  parameters$p_adj <- .ctScoreBH(parameters$p)
  if (q == 1L) parameters$class <- .ctScoreClasses(heldC)
  list(joint = data.frame(covariate = attr(G, 'name'), form = attr(G, 'form'),
    basis = q, LM = lm, df = joint$rank, p = pjoint,
    nmissing = attr(G, 'nmissing'), stringsAsFactors = FALSE),
    parameters = parameters[order(parameters$p, -abs(parameters$estimate /
      parameters$se)), , drop = FALSE])
}

# A covariate that acts above or below a threshold, at a threshold not known
# in advance: the joint test of a step at each candidate cut (the 15th to the
# 85th percentile), and the largest of them, calibrated by one multiplier
# bootstrap over subjects for every cut at once (Andrews, 1993; the split a
# score-guided tree makes, Arnold et al., 2021, 2024). The effects reported
# are those of being above the best cut, held as for a linear covariate.
.ctScoreThreshold <- function(pieces, items, z, name, variance = 'robust',
  nboot = 1999L, ncuts = 15L) {
  if (is.factor(z) || is.character(z) || is.logical(z)) stop("Covariate '", name,
    "' is not numeric; a threshold needs an ordered covariate.", call. = FALSE)
  z <- as.numeric(z)
  ok <- !is.na(z)
  cuts <- unique(stats::quantile(z[ok], seq(0.15, 0.85, length.out = ncuts),
    names = FALSE, type = 1))
  cuts <- cuts[vapply(cuts, function(c) {
    a <- sum(z[ok] > c); a >= 5 && sum(ok) - a >= 5 }, logical(1))]
  if (!length(cuts)) stop("Covariate '", name, "' has too few distinct values ",
    "for a threshold.", call. = FALSE)
  fits <- lapply(cuts, function(c) {
    G <- matrix(as.numeric(z > c), ncol = 1L, dimnames = list(NULL, 'step'))
    G[ok, ] <- G[ok, ] - mean(G[ok, ]); G[!ok, ] <- 0
    attr(G, 'name') <- name; attr(G, 'form') <- 'threshold'
    attr(G, 'nmissing') <- sum(!ok)
    list(G = G)
  })
  # The per-subject projected contributions of each cut, for the sup.
  contrib <- lapply(fits, function(f) .ctScoreCovariateContributions(pieces, items, f$G))
  stat <- function(C) {
    v <- .ctScorePinv(.ctScoreSpread(C))
    u <- colSums(C)
    as.numeric(t(u) %*% v$inv %*% u)
  }
  lms <- vapply(contrib, stat, 0)
  best <- which.max(lms)
  # The best cut's effects and their tests, with the bootstrap the scan skipped.
  out <- .ctScoreCovariate(pieces, items, fits[[best]]$G, variance, nboot)
  out$joint$LM <- lms[best]
  out$joint$p <- .ctScoreBoot(contrib, nboot)
  out$joint$cut <- cuts[best]
  out
}

# Each subject's projected contribution to a covariate's effects: the
# columns whose sum is the joint test's U*.
.ctScoreCovariateContributions <- function(pieces, items, G) {
  S <- pieces$S
  q <- ncol(G); k <- length(items)
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  C <- do.call(cbind, lapply(items, function(it) G * it$s))
  Idt <- do.call(rbind, lapply(items, function(it) -crossprod(G, it$H)))
  C - S %*% t(Idt %*% Iinv)
}

# Relaxations the data cannot tell apart at first order: those whose projected
# per-subject contributions are nearly collinear (|correlation| >= 0.9) would
# each explain the same signal, as two models that are equivalent locally
# (MacCallum et al., 1993). Connected groups get a letter; a relaxation alike
# to no other gets NA.
.ctScoreClasses <- function(C, threshold = 0.9) {
  C <- as.matrix(C)
  k <- ncol(C)
  out <- rep(NA_character_, k)
  if (k < 2L) return(out)
  sdv <- apply(C, 2, stats::sd)
  r <- suppressWarnings(stats::cor(C))
  r[!is.finite(r)] <- 0
  adj <- abs(r) >= threshold & outer(sdv > 0, sdv > 0)
  diag(adj) <- FALSE
  group <- rep(0L, k); g <- 0L
  for (a in seq_len(k)) {
    if (group[a] || !any(adj[a, ])) next
    g <- g + 1L
    todo <- a
    while (length(todo)) {
      b <- todo[1]; todo <- todo[-1]
      if (group[b]) next
      group[b] <- g
      todo <- c(todo, which(adj[b, ] & !group))
    }
  }
  out[group > 0] <- LETTERS[((group[group > 0] - 1L) %% 26L) + 1L]
  out
}

# Whether parameters are stable across the data, ordered by time and by each
# latent's predicted state: Nyblom's (1989) test from each data row's score.
# Under the model a subject's row scores are a martingale difference
# sequence, so their partial sums in any order fixed in advance -- by time, or
# by the state predicted from earlier rows -- behave as a Brownian motion; at
# the estimate they are pulled to zero at the end, a bridge. A parameter that
# drifts over the study, or differs between low and high levels of a state,
# bends its partial sums away from the bridge. L_j is the mean squared partial
# sum of parameter j's score over its total, the Cramer-von Mises form.
#
# The estimation effect is removed by projecting the partial sums with the
# information accumulated up to each row, H_t H_N^-1 S_N, with H the rows'
# outer products -- not with the share of rows, which is right only when the
# information accrues evenly along the ordering, and along a state it does
# not (a drift is informed by large states). That needs every parameter's
# row scores; above `maxall` parameters only the tested ones are used, and the
# others' estimation is ignored. Calibration: a wild bootstrap with a normal
# weight per row, which keeps the rows' variance and nothing else -- what the
# null has, since the scores are uncorrelated over rows (Hansen 1992; Zeileis
# & Hornik 2007 for the fluctuation-test family).
.ctScoreInstability <- function(fit, est, ptab, parameters = NULL, nboot = 1999L,
  cores = 2L, maxall = 48L) {
  tested <- which(ptab$kind == 'model' & !ptab$matrix %in% c('T0MEANS', 'T0VAR'))
  if (!is.null(parameters)) tested <- intersect(tested, which(ptab$name %in% parameters))
  if (!length(tested)) stop("No parameters to test for stability.", call. = FALSE)
  p <- length(est)
  all <- if (p <= maxall) seq_len(p) else tested
  R <- .ctBackendRowScores(fit, est, diag(p)[, all, drop = FALSE], cores)
  err <- .ctAdequacyErrors(fit)
  okrow <- apply(is.finite(R), 1, all)
  R[!okrow, ] <- 0
  Tc <- match(tested, all)
  k <- length(Tc)
  N <- sum(okrow)
  latent <- .ctFitModelObject(fit)$latentNames
  nl <- min(ncol(err$eta), length(latent))
  keys <- c(list(time = err$time), stats::setNames(lapply(seq_len(nl), function(j)
    err$eta[, j]), paste0('state ', latent[seq_len(nl)])))
  rows <- list(); joint <- list()
  for (o in names(keys)) {
    # Ties broken at random: the data's own order within a tie is the
    # subjects', which a statistic over the ordering would read as structure.
    idx <- order(keys[[o]], stats::runif(length(keys[[o]])))
    idx <- idx[okrow[idx]]
    s <- R[idx, , drop = FALSE]
    H <- crossprod(s)
    Hinv <- .ctScorePinv(H)$inv
    VT <- H[Tc, Tc, drop = FALSE]
    VTinv <- .ctScorePinv(VT)$inv
    vd <- diag(VT)
    sT <- s[, Tc, drop = FALSE]
    bridge <- function(sm) {
      B <- Hinv %*% colSums(sm)
      matrix(apply(sm[, Tc, drop = FALSE] - sT * as.numeric(s %*% B), 2, cumsum), ncol = k)
    }
    statistic <- function(St) c(colSums(St^2) / (N * vd), sum((St %*% VTinv) * St) / N)
    observed <- statistic(bridge(s))
    draws <- vapply(seq_len(nboot), function(b) statistic(bridge(s * stats::rnorm(nrow(s)))),
      observed)
    draws <- matrix(draws, nrow = length(observed))
    pv <- (1 + rowSums(draws >= observed)) / (1 + nboot)
    one <- data.frame(ordering = o, parameter = ptab$name[tested], L = observed[seq_len(k)],
      p = pv[seq_len(k)], stringsAsFactors = FALSE)
    one$p_adj <- .ctScoreBH(one$p)
    rows[[o]] <- one
    joint[[o]] <- data.frame(ordering = o, L = observed[k + 1L], df = k, p = pv[k + 1L],
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  list(parameters = out[order(out$p, -out$L), , drop = FALSE],
    joint = do.call(rbind, joint),
    note = if (p > maxall) paste0("Only the ", k, " tested parameters' row scores were ",
      "used; the other parameters' estimation is ignored, which can move the ",
      "level along a state.") else NULL)
}

# The fitted model with a term `d * eta_k * eta_m` in each latent equation's
# CINT, for k <= m (and degree three with `order = 3`), every d a new PARS
# parameter with an identity transform.
.ctScoreStateModel <- function(model, order = 2L) {
  latent <- model$latentNames
  nl <- length(latent)
  # The monomials, each as its sorted state indices: the first is the state
  # the drift coefficient multiplies, the rest what it varies with.
  monomials <- list()
  for (k in seq_len(nl)) for (m in k:nl) {
    monomials[[length(monomials) + 1L]] <- c(k, m)
    if (order >= 3L) for (l in m:nl) monomials[[length(monomials) + 1L]] <- c(k, m, l)
  }
  terms <- list()
  for (i in seq_len(nl)) for (mono in monomials) {
    terms[[length(terms) + 1L]] <- data.frame(name = paste0('ctscoreterm',
      length(terms) + 1L), equation = latent[i], i = i, k = mono[1L],
      by = paste(latent[mono[-1L]], collapse = '*'),
      term = paste(latent[mono], collapse = '*'),
      other = if (length(mono) == 2L && mono[1L] != mono[2L]) mono[2L] else NA_integer_,
      stringsAsFactors = FALSE)
  }
  terms <- do.call(rbind, terms)
  model <- .ctAddCintTerms(model, terms$i, terms$name,
    lapply(strsplit(terms$term, '*', fixed = TRUE), match, latent))
  list(model = model, terms = terms)
}

# Add terms `name * eta_k * eta_m ...` to CINT[i], each `name` a new PARS
# parameter (identity transform unless given), for the score tests at zero and
# for a refit that keeps one. A CINT cell that is a free parameter moves into
# PARS -- transform, variation and predictor effects and all -- so its raw
# coordinate is unchanged; a fixed cell keeps its value; an expression keeps its
# form. The engine's Jacobian is derived from DRIFT * state + CINT
# (`ctJacobian`), so a term here is the same model as the drift coefficient
# varying with the state.
.ctAddCintTerms <- function(model, i, names, monomials, transform = 'param') {
  pars <- model$pars
  latent <- model$latentNames
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  nextpars <- max(c(0L, as.integer(pars$row[pars$matrix %in% 'PARS'])))
  template <- pars[which(pars$matrix %in% 'CINT')[1L], , drop = FALSE]
  parsrow <- function(name, from = NULL) {
    r <- if (is.null(from)) template else from
    r$matrix <- 'PARS'; r$row <- nextpars <<- nextpars + 1L; r$col <- 1L
    r$param <- name
    if (is.null(from)) {
      r$value <- NA; r$transform <- transform; r$indvarying <- FALSE
      if ('sdscale' %in% names(r)) r$sdscale <- 1
      for (e in effects) r[[e]] <- 'FALSE'
      for (l in levelflags) r[[l]] <- FALSE
    }
    r
  }
  added <- list()
  for (eq in unique(i)) {
    cell <- which(pars$matrix %in% 'CINT' & pars$row == eq)
    if (length(cell) != 1L) stop("The model has no CINT cell for latent ", eq,
      ".", call. = FALSE)
    these <- which(i == eq)
    sum_i <- paste(vapply(these, function(t) paste(c(names[t], latent[monomials[[t]]]),
      collapse = ' * '), character(1)), collapse = ' + ')
    for (t in these) added[[length(added) + 1L]] <- parsrow(names[t])
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
  model
}

# The fixed cells a model could free, as terms at zero: DRIFT cross effects,
# DIFFUSION and MANIFESTVAR covariances, LAMBDA cross loadings and
# TDPREDEFFECT cells fixed at zero. Each becomes a parameter with ctsem's
# default transform for its cell, kept only where that transform puts a raw
# zero at zero, so the extended model is the fit at zero. CINT and
# MANIFESTMEANS are left out: with the other free, either is identified only
# by the transient from the initial state.
.ctScoreFixedModel <- function(model) {
  pars <- model$pars
  latent <- model$latentNames; manifest <- model$manifestNames
  tdpred <- model$TDpredNames
  ct <- isTRUE(model$continuoustime)
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  value <- suppressWarnings(as.numeric(pars$value))
  free <- is.na(pars$param) | !nzchar(pars$param)
  candidate <- which(free & !is.na(value) & value == 0 & (
    (pars$matrix %in% 'DRIFT' & pars$row != pars$col) |
      (pars$matrix %in% c('DIFFUSION', 'MANIFESTVAR') & pars$row > pars$col) |
      pars$matrix %in% c('LAMBDA', 'TDPREDEFFECT')))
  ns <- asNamespace('ctsem')
  dims <- function(m) switch(m, DRIFT = list(latent, latent),
    DIFFUSION = list(latent, latent), MANIFESTVAR = list(manifest, manifest),
    LAMBDA = list(manifest, latent), TDPREDEFFECT = list(latent, tdpred))
  terms <- list()
  for (k in candidate) {
    m <- pars$matrix[k]; r <- pars$row[k]; c <- pars$col[k]
    tf <- .ctModelDefaultFreePar(m, r, c, ct)$transform
    atzero <- try(eval(parse(text = tf), list(param = 0), ns), silent = TRUE)
    if (inherits(atzero, 'try-error') || !isTRUE(abs(atzero) < 1e-12)) next
    d <- dims(m)
    name <- paste0('ctscorefixed', length(terms) + 1L)
    pars$param[k] <- name
    pars$value[k] <- NA
    pars$transform[k] <- tf
    pars$indvarying[k] <- FALSE
    for (e in effects) pars[[e]][k] <- 'FALSE'
    for (l in levelflags) pars[[l]][k] <- FALSE
    terms[[length(terms) + 1L]] <- data.frame(name = name, matrix = m, row = r,
      col = c, cell = paste0(m, '[', d[[1]][r], ',', d[[2]][c], ']'),
      suggested = paste0(tolower(m), '_', d[[1]][r], '_', d[[2]][c]),
      transform = tf, stringsAsFactors = FALSE)
  }
  model$pars <- pars
  list(model = model, terms = if (length(terms)) do.call(rbind, terms) else NULL)
}

# The fitted model extended at zero by every term of the families asked for --
# drift terms written into CINT, fixed cells freed -- prepared as the fit was,
# evaluated with the fitted values and the terms at zero, which must reproduce
# the fit; then its gradient, its Hessian columns along the terms, and for the
# robust variance each unit's gradient, whose term columns projected on the
# fitted parameters are the contributions. One model for both families, so
# one compile of the extended shape.
.ctScoreExtended <- function(fit, est, families, information = NULL, order = 2L,
  variance = 'robust', nboot = 1999L, cores = 2L) {
  spec <- fit$model_spec
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) stop("This fit does not carry the model it was built from, ",
    "so the extended model cannot be prepared. Refit with this version of ctsem.",
    call. = FALSE)
  model <- base
  state <- fixed <- NULL
  if ('statedependence' %in% families) {
    state <- .ctScoreStateModel(model, order)
    model <- state$model
  }
  if ('fixed' %in% families) {
    fixed <- .ctScoreFixedModel(model)
    model <- fixed$model
  }
  names_new <- c(state$terms$name, fixed$terms$name)
  if (!length(names_new)) return(list(statedependence = NULL, fixed = NULL))
  args <- .ctFitReplayArgs(fit)
  args$model <- model
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
    "State dependence and fixed cells cannot be tested for this model.", call. = FALSE)
  if (!setequal(setdiff(nm2, nm1), names_new)) stop("The extended model added ",
    "parameters other than its terms: ", paste(setdiff(setdiff(nm2, nm1), names_new),
      collapse = ', '), call. = FALSE)
  x2 <- setNames(numeric(length(nm2)), nm2)
  x2[nm1] <- est
  value1 <- ctJuliaEvaluate(fit, est, gradient = FALSE, cores = cores)$value
  ev <- ctJuliaEvaluate(spec2, x2, gradient = TRUE, cores = cores)
  # With priors = TRUE every raw parameter has a normal(0, 1) prior, the new
  # terms too: at zero each adds log dnorm(0) to the objective, nothing to its
  # gradient and -1 to its own curvature. The terms are tested on the
  # likelihood, so that offset is expected and the curvature is taken back out
  # below. Before this, a four-latent fit with priors refused at 40 terms x
  # -0.919 = -36.76, which read as a model the rebuild had got wrong.
  priored <- isTRUE(fit$args$resolved$priors)
  offset <- if (priored) length(names_new) * stats::dnorm(0, log = TRUE) else 0
  if (!isTRUE(abs(ev$value - value1 - offset) <= 1e-6 * max(1, abs(value1)))) stop(
    "The extended model does not reproduce the fit at zero (log likelihood ",
    signif(ev$value - offset, 8), " against ", signif(value1, 8), "), so its scores are ",
    "not the fit's. State dependence and fixed cells cannot be tested for this model.",
    call. = FALSE)
  grad <- setNames(as.numeric(ev$gradient), nm2)
  theta <- match(nm1, nm2)
  delta <- match(names_new, nm2)
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(spec2, hybrid = TRUE)
  Hc <- .ctBackendWithMaxChunks(cores, as.matrix(.ctBackendJuliaValue(
    module$ctsem_hessian_columns(objective, .ctJuliaNumericVector(as.numeric(x2)),
      .ctJuliaVector(as.integer(delta))))))
  if (priored) Hc[cbind(delta, seq_along(delta))] <- Hc[cbind(delta, seq_along(delta))] + 1
  if (is.null(information)) information <- -.ctBackendHessian(fit, est)
  Iinv <- .ctScorePinv(information)$inv
  Idt <- -t(Hc[theta, , drop = FALSE])
  Idd <- -Hc[delta, , drop = FALSE]
  Idd <- (Idd + t(Idd)) / 2
  A <- Idt %*% Iinv
  Ustar <- grad[delta] - as.numeric(A %*% grad[theta])
  Vmodel <- Idd - A %*% t(Idt)
  Vmodel <- (Vmodel + t(Vmodel)) / 2
  V <- Vmodel
  Cstar <- NULL
  if (variance == 'robust') {
    units <- .ctBackendWithMaxChunks(cores, .ctJuliaGet(module$ctsem_subject_gradients(
      objective, .ctJuliaVector(as.numeric(x2)))))
    units <- as.matrix(units$scores)
    Cstar <- units[, delta, drop = FALSE] - units[, theta, drop = FALSE] %*% t(A)
    V <- .ctScoreSpread(Cstar)
  }
  v <- diag(V)
  vm <- diag(Vmodel)
  per <- data.frame(LM = ifelse(v > 0, Ustar^2 / v, NA_real_),
    # One term at a time, the fitted parameters free: the expected parameter
    # change. NA where the fitted model's curvature cannot size it -- a strong
    # departure can leave the projected information indefinite.
    estimate = ifelse(vm > 0, Ustar / vm, NA_real_),
    se = ifelse(vm > 0, sqrt(pmax(v, 0)) / vm, NA_real_))
  per$p <- stats::pchisq(per$LM, 1, lower.tail = FALSE)
  jointtest <- function(idx, label) {
    j <- .ctScorePinv(V[idx, idx, drop = FALSE])
    stat <- as.numeric(t(Ustar[idx]) %*% j$inv %*% Ustar[idx])
    p <- if (!j$rank) NA_real_ else if (!is.null(Cstar))
      .ctScoreBoot(Cstar[, idx, drop = FALSE], nboot) else
      stats::pchisq(stat, j$rank, lower.tail = FALSE)
    data.frame(set = label, LM = stat, df = j$rank, p = p, stringsAsFactors = FALSE)
  }
  out <- list(statedependence = NULL, fixed = NULL)
  latent <- base$latentNames
  if (!is.null(state$terms)) {
    idx <- seq_len(nrow(state$terms))
    terms <- state$terms
    table <- data.frame(term = terms$term, equation = terms$equation,
      reads = paste0('DRIFT[', latent[terms$i], ',', latent[terms$k], '] varies with ',
        terms$by, ifelse(is.na(terms$other), '', paste0(', or DRIFT[',
          latent[terms$i], ',', latent[pmax(1L, terms$other)], '] with ',
          latent[terms$k]))),
      per[idx, c('LM', 'p', 'estimate', 'se')], stringsAsFactors = FALSE)
    table$p_adj <- .ctScoreBH(table$p)
    if (!is.null(Cstar)) table$class <- .ctScoreClasses(Cstar[, idx, drop = FALSE])
    table$i <- terms$i; table$k <- terms$k
    table$monomial <- vapply(seq_len(nrow(terms)), function(r)
      paste(match(strsplit(terms$term[r], '*', fixed = TRUE)[[1]], latent), collapse = ','),
      character(1))
    joint <- do.call(rbind, c(lapply(unique(terms$equation), function(e)
      jointtest(idx[terms$equation == e], e)), list(jointtest(idx, 'all'))))
    names(joint)[1] <- 'equation'
    out$statedependence <- list(terms = table[order(table$p, -table$LM), , drop = FALSE],
      joint = joint)
  }
  if (!is.null(fixed$terms)) {
    idx <- length(state$terms$name) + seq_len(nrow(fixed$terms))
    terms <- fixed$terms
    table <- data.frame(cell = terms$cell, matrix = terms$matrix, row = terms$row,
      col = terms$col, per[idx, c('LM', 'p', 'estimate', 'se')],
      suggested = terms$suggested, transform = terms$transform, stringsAsFactors = FALSE)
    table$p_adj <- .ctScoreBH(table$p)
    if (!is.null(Cstar)) table$class <- .ctScoreClasses(Cstar[, idx, drop = FALSE])
    out$fixed <- list(cells = table[order(table$p, -table$LM), , drop = FALSE],
      joint = jointtest(idx, 'all'))
  }
  out
}

#' @export
print.ctScoreTest <- function(x, digits = 3, top = 10, ...) {
  # A bootstrap p-value cannot go below 1 / (nboot + 1); at that floor it
  # says only that, and the table is ordered by the statistic instead.
  floor <- if (identical(x$info$variance, 'robust') && !is.null(x$info$nboot))
    1 / (x$info$nboot + 1) else 0
  fmt <- function(d) {
    d <- utils::head(d, top)
    atfloor <- if ('p' %in% names(d)) !is.na(d$p) & d$p <= floor * (1 + 1e-9) else NULL
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], function(v) signif(v, digits))
    if (any(atfloor)) {
      d$p <- format(d$p)
      d$p[atfloor] <- paste0('<', signif(floor, 2))
    }
    print(d, row.names = FALSE)
  }
  cat("Score tests at the fitted model (", x$info$route, " route, ",
    x$info$nsubjects, " subjects)\n", sep = "")
  for (e in names(x$info$errors)) cat("Not tested, ", e, ": ", x$info$errors[[e]], "\n", sep = "")
  withclass <- function(d, cols) if (any(!is.na(d$class))) c(cols, 'class') else cols
  if (!is.null(x$heterogeneity)) {
    cat("\nUnmodelled heterogeneity (one-sided; sd is the one-step estimate, raw scale):\n")
    fmt(x$heterogeneity[, withclass(x$heterogeneity, c('parameter', 'z', 'p', 'p_adj',
      'sd', 'sd_natural'))])
  }
  if (!is.null(x$dimension)) {
    d <- x$dimension
    cat("\nA new dimension of variation: largest eigenvalue ", signif(d$lambda, digits),
      ", p = ", signif(d$p, digits), "; direction ",
      paste0(names(d$direction)[seq_len(min(5, length(d$direction)))], ' ',
        signif(d$direction[seq_len(min(5, length(d$direction)))], 2), collapse = ', '),
      "\n", sep = "")
  }
  for (c in x$covariates) {
    j <- c$joint
    cat("\nCovariate ", j$covariate, " (", j$form,
      if (!is.null(j$cut)) paste0(", above ", signif(j$cut, digits)), "): joint LM = ",
      signif(j$LM, digits), ", df = ", j$df, ", p = ", signif(j$p, digits), "\n", sep = "")
    fmt(c$parameters[, withclass(c$parameters, c('parameter', 'estimate', 'se',
      'estimate_natural', 'p', 'p_adj'))])
  }
  if (!is.null(x$statedependence)) {
    cat("\nState dependence of the drift:\n")
    fmt(x$statedependence$joint)
    cat("\n")
    fmt(x$statedependence$terms[, withclass(x$statedependence$terms, c('term',
      'equation', 'LM', 'p', 'p_adj', 'estimate', 'se', 'reads'))])
  }
  if (!is.null(x$instability)) {
    cat("\nStability across the data (Nyblom), by ordering:\n")
    fmt(x$instability$joint)
    cat("\n")
    fmt(x$instability$parameters[, c('ordering', 'parameter', 'L', 'p', 'p_adj')])
    if (!is.null(x$instability$note)) cat(x$instability$note, "\n")
  }
  if (!is.null(x$fixed)) {
    j <- x$fixed$joint
    cat("\nFixed cells, freed: joint LM = ", signif(j$LM, digits), ", df = ", j$df,
      ", p = ", signif(j$p, digits), "\n", sep = "")
    fmt(x$fixed$cells[, withclass(x$fixed$cells, c('cell', 'LM', 'p', 'p_adj',
      'estimate', 'se'))])
  }
  invisible(x)
}
