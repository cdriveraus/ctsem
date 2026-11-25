#' Drift diffusion model for choice and response time
#'
#' Builds a ctsem model of choices and response times by the drift diffusion
#' model, with the drift as a latent process. The model is an ordinary
#' \code{\link{ctModel}} with one choice and response time indicator
#' (\code{manifesttype = 5}), so it is fitted with \code{\link{ctFit}} and can
#' be changed or extended like any other model, for example with other
#' indicators of the same latent process.
#'
#' @param name Name of the data column holding each trial's response time in
#'   seconds, negative for a response at the lower boundary. With columns
#'   \code{rt} and \code{correct}, \code{dat$rt_signed <- ifelse(dat$correct ==
#'   1, dat$rt, -dat$rt)} makes it for accuracy coding.
#' @param drift,boundary,ndt One sided formulas of trial level predictors of
#'   the drift, the log boundary separation and the log nondecision time, such
#'   as \code{~ speed + difficulty}. \code{~ 1} for none. Predictors are numeric
#'   columns of the data, entered as time dependent predictors; code a factor
#'   as dummy columns. The drift's predictors shift the drift on their trial
#'   only.
#' @param bias The relative start point, as a proportion of the boundary
#'   separation from the lower boundary: a number to fix it, or a one sided
#'   formula to estimate it, with predictors on the logit scale.
#' @param random Which of \code{'drift'}, \code{'boundary'}, \code{'ndt'} and
#'   \code{'bias'} vary over subjects.
#' @param dynamic Whether the drift changes over the session. \code{FALSE}
#'   keeps each subject's drift constant. \code{TRUE} makes it a mean
#'   reverting (Ornstein-Uhlenbeck) process in the data's time column around
#'   the subject's level, starting from its stationary distribution. The time
#'   column sets the scale of its rate: trial number, or seconds.
#' @param sv Whether to estimate trial to trial variability of the drift, as a
#'   standard deviation.
#' @param ... Passed to \code{\link{ctModel}}, such as \code{id}, \code{time}
#'   or \code{TIpredNames}.
#'
#' @return A \code{ctStanModel} object.
#' @details The drift's parameters are named \code{drift_mean} (its mean, or
#'   with \code{dynamic = TRUE} \code{drift_cint}, \code{drift_auto} and
#'   \code{drift_diffusion} for the process's intercept, auto effect and
#'   diffusion), \code{boundary}, \code{ndt} and \code{bias}, and each
#'   predictor's effect \code{<parameter>_<predictor>}. The diffusion
#'   coefficient of the decision process is 1.
#' @seealso \code{\link{ctModel}}, whose \code{WIENER} argument sets the
#'   boundary, nondecision time and bias of a choice and response time
#'   indicator directly.
#' @examples
#' m <- ctModelDDM(name = 'rt_signed', drift = ~ difficulty, boundary = ~ speed,
#'   dynamic = TRUE)
#' print(m)
#' @export
ctModelDDM <- function(name = 'rt', drift = ~ 1, boundary = ~ 1, ndt = ~ 1,
    bias = 0.5, random = c('drift', 'boundary', 'ndt'), dynamic = FALSE,
    sv = FALSE, ...) {
  random <- match.arg(random, c('drift', 'boundary', 'ndt', 'bias'),
    several.ok = TRUE)
  biasformula <- inherits(bias, 'formula')
  if (!biasformula && !(is.numeric(bias) && length(bias) == 1L && bias > 0 &&
      bias < 1)) stop('bias is a number between 0 and 1, or a one sided formula',
    call. = FALSE)
  if (!biasformula && 'bias' %in% random) stop("a fixed bias cannot vary over ",
    "subjects; give bias = ~ 1 to estimate it", call. = FALSE)
  preds <- list(drift = .ctDDMPredictors(drift, 'drift'),
    boundary = .ctDDMPredictors(boundary, 'boundary'),
    ndt = .ctDDMPredictors(ndt, 'ndt'),
    bias = if (biasformula) .ctDDMPredictors(bias, 'bias') else character())
  shift <- function(par) {
    x <- preds[[par]]
    if (!length(x)) return('')
    paste(paste0(par, '_', x, ' * ', x), collapse = ' + ')
  }
  vary <- function(par) if (par %in% random) 'TRUE' else 'FALSE'
  pars <- character()
  cell <- function(par, transform) {
    if (!length(preds[[par]])) return(paste0(par, '|', transform, '|', vary(par)))
    pars <<- c(pars, paste0(par, '|', transform, '|', vary(par)))
    switch(par,
      bias = paste0('inv_logit(', par, ' + ', shift(par), ')'),
      paste0(par, ' * exp(', shift(par), ')'))
  }
  # A bias with predictors is on the logit scale in its cell; without them
  # the WIENER default transform is the inverse logit.
  wiener <- c(cell('boundary', 'exp(param)'),
    cell('ndt', paste0('exp(param + ', signif(log(0.1), 15), ')')),
    if (!biasformula) as.character(bias) else if (length(preds$bias))
      cell('bias', 'param') else paste0('bias||', vary('bias')))
  effects <- unlist(lapply(names(preds), function(p) if (length(preds[[p]]))
    paste0(p, '_', preds[[p]], '||FALSE')))
  td <- unique(unlist(preds))
  args <- list(type = 'ct', manifestNames = name, manifesttype = 5,
    latentNames = 'drift', LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(if (length(preds$drift)) shift('drift') else 0),
    MANIFESTVAR = matrix(if (sv) 'sv|log1p_exp(param)|FALSE' else 0),
    WIENER = wiener)
  if (length(td)) {
    args$TDpredNames <- td
    args$TDPREDEFFECT <- matrix(0, 1, length(td))
  }
  if (length(c(pars, effects))) args$PARS <- c(pars, effects)
  if (dynamic) {
    args$DRIFT <- matrix('drift_auto')
    args$DIFFUSION <- matrix('drift_diffusion')
    args$CINT <- matrix(paste0('drift_cint||', vary('drift')))
    args$T0MEANS <- matrix('stationary')
    args$T0VAR <- matrix('stationary')
  } else {
    args$DRIFT <- matrix(0)
    args$DIFFUSION <- matrix(0)
    args$CINT <- matrix(0)
    args$T0MEANS <- matrix(paste0('drift_mean||', vary('drift')))
    args$T0VAR <- matrix(1e-3)
  }
  extra <- list(...)
  clash <- intersect(names(extra), names(args))
  if (length(clash)) stop(paste(clash, collapse = ', '), ' is set by ',
    'ctModelDDM(); build the model with ctModel() to set it yourself',
    call. = FALSE)
  suppressMessages(do.call(ctModel, c(args, extra)))
}

# The predictor names of a one sided formula, main effects only.
.ctDDMPredictors <- function(f, what) {
  if (!inherits(f, 'formula') || length(f) != 2L) stop(what, ' must be a one ',
    'sided formula, such as ~ 1 or ~ speed', call. = FALSE)
  labels <- attr(stats::terms(f), 'term.labels')
  if (any(grepl(':', labels, fixed = TRUE))) stop('interactions in ', what,
    ' are not supported; make a column of the product', call. = FALSE)
  if (any(labels != make.names(labels))) stop(what, ' takes column names ',
    'only, not expressions: ', paste(labels[labels != make.names(labels)],
      collapse = ', '), call. = FALSE)
  labels
}
