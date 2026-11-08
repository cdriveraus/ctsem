#' Likelihood-ratio tests between nested ctsem fits
#'
#' \code{anova()} on two or more optimised fits of nested models to the same
#' data gives each fit's number of parameters, AIC, BIC, log likelihood and -2
#' log likelihood, and tests each fit against the one above it: the chi-square
#' statistic (twice the gain in log likelihood), its degrees of freedom (the
#' extra parameters) and its p-value. The fits may come from either backend.
#'
#' The fits are ordered by their number of parameters, and each is tested
#' against the one with the next fewer, so the models must be nested in that
#' order; \code{anova()} cannot check that. A negative statistic, which warns,
#' means a model fits worse than one nested in it: the models are not nested,
#' or an optimisation stopped short.
#'
#' The p-value is the chi-square distribution's on the difference in
#' parameters. Where the smaller model fixes a variance at zero, such as a
#' random effect's, the null lies on the boundary of the parameter space and
#' that p-value is conservative (Stram & Lee, 1994); \code{\link{ctScoreTest}}
#' tests a new random effect against the corrected reference without
#' refitting.
#'
#' The log likelihoods are those \code{logLik()} reports. Where a fit's
#' random-effects check or a refined fit integrated the random effects, that
#' is the integrated value, and the heading says so. A fit with priors
#' maximised the posterior, so its log likelihood is read at the posterior
#' mode rather than at the likelihood's maximum, and the heading says that too.
#' Sampled fits have no single log likelihood; compare them with
#' \code{\link{ctLOO}}.
#'
#' @param object An optimised fit from \code{\link{ctFit}}.
#' @param ... Further optimised fits, of models nested with \code{object} and
#' fitted to the same data.
#'
#' @return A data frame of class \code{anova}, one row per fit in order of
#' parameters, with columns \code{npar}, \code{AIC}, \code{BIC},
#' \code{logLik}, \code{-2LL}, \code{Chisq}, \code{Df} and
#' \code{Pr(>Chisq)}; the last three are \code{NA} on the first row. Rows are
#' named by the arguments.
#'
#' @references Stram, D. O., & Lee, J. W. (1994). Variance components testing
#' in the longitudinal mixed effects model. \emph{Biometrics}, 50(4),
#' 1171-1177.
#'
#' @seealso \code{\link{ctChisqTest}} for the p-value alone.
#' @name anova.ctJuliaFit
#' @aliases anova.ctStanFit
#' @examples
#' \donttest{
#'     df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)
#'
#'     m1 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0)
#'     m2 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)
#'
#'     if (isTRUE(ctJuliaStatus()$available)) {
#'       f1 <- ctFit(df,m1,cores=1)
#'       f2 <- ctFit(df,m2,cores=1)
#'
#'       test <- anova(f1, f2)
#'       test
#'       test$Chisq[2]
#'     }
#' }
#' @export
anova.ctJuliaFit <- function(object, ...) {
  # Every fit named, as in anova(RI = fitRI, FI = fitFI), leaves `object`
  # missing and all of them in the dots.
  fits <- if (missing(object)) list(...) else list(object, ...)
  .ctAnova(fits, .ctAnovaNames(match.call(), length(fits)))
}

#' @rdname anova.ctJuliaFit
#' @export
anova.ctStanFit <- function(object, ...) {
  fits <- if (missing(object)) list(...) else list(object, ...)
  .ctAnova(fits, .ctAnovaNames(match.call(), length(fits)))
}

# The rows' names: each argument as written (`anova(fitRI, fitFI)`), or its
# name where one is given (`anova(RI = fitRI, FI = fitFI)`), else fit1, fit2,
# ... match.call() lists `object` first and then the dots, the order the fits
# are listed in. A fit handed over by value, as do.call() hands it, is never
# deparsed: that would print the whole object.
.ctAnovaNames <- function(call, n) {
  args <- as.list(call)[-1L]
  given <- names(args)
  if (is.null(given)) given <- rep('', length(args))
  lab <- vapply(args, function(e) if (is.name(e) || is.call(e))
    paste(deparse(e, width.cutoff = 500L), collapse = ' ') else '', '')
  use <- nzchar(given) & given != 'object'
  lab[use] <- given[use]
  if (length(lab) != n || !all(nzchar(lab)) || anyDuplicated(lab) ||
      any(nchar(lab) > 40L)) lab <- paste0('fit', seq_len(n))
  unname(lab)
}

.ctFitSampled <- function(fit) {
  if (.ctFitIsJulia(fit)) !is.null(fit$sample) else
    length(fit$stanfit$stanfit@sim) > 0
}

# TRUE, 'randomCorr' or stan's 1: any prior makes the estimate a posterior
# mode, where the log likelihood is not at its maximum.
.ctFitHasPriors <- function(fit) {
  p <- if (.ctFitIsJulia(fit)) .ctJuliaOr(fit$args$resolved$priors,
    fit$args$priors) else fit$standata$priors
  length(p) == 1L && !is.na(p) && !identical(p, FALSE) &&
    !isTRUE(p == 0)
}

# The table behind anova() and ctChisqTest(): every fit's logLik(), ordered by
# parameters, each tested against the row above.
.ctAnova <- function(fits, names) {
  if (length(fits) < 2L) stop("anova() compares nested fits; give it two or ",
    "more.", call. = FALSE)
  isfit <- vapply(fits, function(f) inherits(f, c('ctJuliaFit', 'ctStanFit')), NA)
  if (!all(isfit)) stop("anova() compares ctsem fits; ",
    paste(names[!isfit], collapse = ', '), " is not one.", call. = FALSE)
  # A likelihood-ratio test needs the log likelihood at a point estimate, which
  # a sampled fit does not have: it has a posterior.
  sampled <- vapply(fits, .ctFitSampled, NA)
  if (any(sampled)) stop(paste(names[sampled], collapse = ' and '),
    if (sum(sampled) > 1L) " are" else " is", " sampled (ctFit(optimize=FALSE) ",
    "or ctFitUncertainty(fit, 'sample')): a sampled fit has a posterior, not a ",
    "single log likelihood. Compare optimized fits, or compare sampled julia ",
    "fits with ctLOO(fit, method = 'psis').", call. = FALSE)
  ll <- lapply(fits, stats::logLik)
  nobs <- vapply(ll, function(l) as.numeric(.ctJuliaOr(attr(l, 'nobs'), NA)), 0)
  nsub <- vapply(fits, function(f) as.numeric(.ctJuliaOr(f$standata$nsubjects, NA)), 0)
  differs <- function(x) length(unique(x[!is.na(x)])) > 1L
  if (differs(nobs) || differs(nsub)) stop("anova() compares fits to the same ",
    "data; these were fitted to ", paste0(nobs, ' observations',
      ifelse(is.na(nsub), '', paste0(' of ', nsub, ' subjects')), ' (',
      names, ')', collapse = ', '), ".", call. = FALSE)

  npar <- vapply(ll, function(l) as.numeric(attr(l, 'df')), 0)
  value <- vapply(ll, as.numeric, 0)
  o <- order(npar)
  npar <- npar[o]; value <- value[o]; ll <- ll[o]; fits <- fits[o]
  names <- names[o]
  chisq <- c(NA, 2 * diff(value))
  df <- c(NA, diff(npar))
  p <- ifelse(!is.na(df) & df > 0,
    stats::pchisq(chisq, pmax(df, 1), lower.tail = FALSE), NA_real_)
  tab <- data.frame(npar = npar,
    AIC = vapply(ll, stats::AIC, 0), BIC = vapply(ll, stats::BIC, 0),
    logLik = value, `-2LL` = -2 * value, Chisq = chisq, Df = df,
    `Pr(>Chisq)` = p, row.names = names, check.names = FALSE)

  worse <- which(!is.na(chisq) & chisq < 0 & df > 0)
  if (length(worse)) warning(paste0(names[worse], " has more parameters than ",
    names[worse - 1L], " but a lower log likelihood", collapse = '; '),
    ": the models are not nested in that order, or an optimisation stopped ",
    "short.", call. = FALSE)

  method <- vapply(fits, function(f) if (.ctFitIsJulia(f))
    as.character(.ctJuliaOr(f$estimate$loglik_method, '')) else '', '')
  notes <- c(
    vapply(which(method == 'checked'), function(i) paste0(names[i],
      ": log likelihood integrated by the random-effects check, ",
      .ctLoglikCheckedError(fits[[i]]$estimate), "."), ''),
    vapply(which(method == 'quadrature'), function(i) paste0(names[i],
      ": log likelihood integrated by the refined fit's quadrature."), ''),
    vapply(which(vapply(fits, .ctFitHasPriors, NA)), function(i)
      paste0(names[i], ": fitted with priors, so its log likelihood is at the ",
        "posterior mode, not at the likelihood's maximum."), ''))
  structure(tab, heading = c(paste0("Likelihood-ratio tests, each fit against ",
    "the one above it", if (length(notes)) "" else "\n"),
    if (length(notes)) c(notes, "")),
    class = c('anova', 'data.frame'))
}

#' Chi Square test wrapper for ctsem fit objects.
#'
#' The p-value of a likelihood-ratio test between two nested model fits, from
#' either backend. \code{\link[=anova.ctJuliaFit]{anova}(fit1, fit2)} gives
#' the whole test: each fit's log likelihood, the chi-square statistic, its
#' degrees of freedom and the p-value.
#'
#' @param fit1 One of the fits to be compared.
#' @param fit2 Second fit to be compared. The fit with more parameters is
#' tested against the other, whichever order they come in.
#'
#' @return The p-value: that of the fit with more parameters in
#' \code{anova(fit1, fit2)}. \code{NA} where the two have the same number of
#' parameters.
#' @export
#'
#' @examples
#' \donttest{
#'     df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)
#'
#'     m1 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0)
#'     m2 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)
#'
#'     if (isTRUE(ctJuliaStatus()$available)) {
#'       f1 <- ctFit(df,m1,cores=1)
#'       f2 <- ctFit(df,m2,cores=1)
#'
#'       ctChisqTest(f1,f2)
#'     }
#' }

ctChisqTest<-function(fit1,fit2){
  tab <- .ctAnova(list(fit1, fit2), c('fit1', 'fit2'))
  tab[['Pr(>Chisq)']][2L]
}
