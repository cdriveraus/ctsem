#' ctRawParnames
#'
#' Gets internal stan parameter names of a ctStanFit object sampled via stan based on specified substrings.
#' \code{ctStanParnames} is maintained as a backward-compatible alias.
#'
#' @param x ctStanFit object, sampled (\code{ctFit(backend = 'stan', optimize = FALSE)}).
#' @param substrings vector of character strings, parameter names of the stan model
#' containing any of these strings will be returned. Useful strings may be 'pop_' for
#' population means, 'popsd' for population standard deviations,
#'  or specific combinations such as 'pop_DRIFT' for the population
#' means of temporal dynamics parameters
#' @return vector of character strings.
#' @examples
#' \donttest{
#' sunspots<-sunspot.year
#' sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#' id <- 1
#' time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#'
#' #setup model
#' ssmodel <- ctModel(type='ct', n.latent=2, n.manifest=1,
#'  manifestNames='sunspots',
#'  latentNames=c('ss_level', 'ss_velocity'),
#'  LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
#'  DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
#'  MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
#'  MANIFESTVAR=diag(0,1), #As per original spec
#'  CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'  DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#'
#' #fit
#' ssfit <- ctFit(datalong, ssmodel, optimize=FALSE, priors=TRUE,
#'   sampleControl=list(iter=2, chains=1), backend='stan')
#' ctRawParnames(ssfit,substrings=c('pop_','popsd'))
#' }
#'
#' @aliases ctStanParnames
#' @export
ctRawParnames <- function(x,substrings=c('pop_','popsd')){
  if(length(x$stanfit$stanfit@sim)==0) stop('Doesnt contain sampled stanfit object')
  out<-c()
  for(subsi in substrings){
    out<- c(out, x$stanfit$stanfit@model_pars[grep(paste0('^',subsi),x$stanfit$stanfit@model_pars)])
  }
  return(out)
}

#' @export
ctStanParnames <- ctRawParnames



#'ctDiscretePars
#'
#'Calculate model implied regressions for a sequence of time intervals (if ct) or steps (if dt) based on
#'a ctStanFit object, for specified subjects. Wrap with print() when used inside for loops!
#'\code{ctStanDiscretePars} is maintained as a backward-compatible alias.
#'
#'@param fit model fit from \code{\link{ctFit}}
#'@param ctstanfitobj Deprecated. Use \code{fit}.
#'@param subjects Either 'popmean', to use the population mean parameter, or the
#'subjects to use, as ids in the data (see \code{realid}).
#'@param realid If TRUE (the default), subjects are read as the ids in the data,
#' and integers that are not ids there are read as ctsem's numbering of the
#' subjects, 1 to N, with a message. If FALSE, they are read as that numbering.
#'@param times Numeric vector of positive values, discrete time parameters will be calculated for each. If the fit
#'object is a discrete time model, these should be positive integers.
#'@param nsamples Number of samples from the stanfit to use for plotting. Higher values will
#'increase smoothness / accuracy, at cost of plotting speed. Values greater than the total
#'number of samples will be set to total samples.
#'@param impulseType What a "one unit change in process c" is taken to bring
#'with it. Every variant is \code{dtDRIFT(t) \%*\% C} for a different companion
#'matrix \code{C}, whose column c says what else moves.
#'\describe{
#'  \item{\code{'unit'} (the default)}{\code{C = I}. An impulse of one unit
#'    (one standard deviation when \code{standardise=TRUE}) is applied to
#'    process c alone and nothing else moves, because the change was imposed.
#'    This is the \emph{partial} regression \code{E[x(t+u)|x(t)]}, and the
#'    only variant that is a property of the dynamics alone.}
#'  \item{\code{'observed'}}{\code{C = Sigma diag(Sigma)^-1}
#'    with \code{Sigma = asymDIFFUSIONcov}. Process c is \emph{observed} one
#'    unit above expectation and the others are not held fixed, so they move by
#'    \code{Sigma_rc/Sigma_cc}. This is the \emph{simple} regression -- what
#'    regressing the observed data on x_c alone recovers.}
#'  \item{\code{'shock'}}{\code{C = Q diag(Q)^-1} with
#'    \code{Q = DIFFUSIONcov}. One hypothetical system noise innovation arrives
#'    in c, with companion innovations \code{Q_rc/Q_cc}. An impulse response
#'    rather than a regression.}
#'  \item{\code{'orthogonal'}}{From the Cholesky factor of DIFFUSION, so
#'    shocks are uncorrelated. Depends on the order of the latent processes,
#'    which is a modelling assumption rather than a property of the fit.}
#'}
#'The companion matrices are asymmetric, and must be: a one unit move in a
#'wide-spread process implies more movement in a narrow one than the reverse.
#'\code{'observed'} and \code{'shock'} use different covariance matrices
#'because they ask different questions, and coincide only under isotropic decay
#'with no cross effects.
#'@param observational Deprecated. Use \code{impulseType}: \code{FALSE} is
#'\code{'unit'} and \code{TRUE} is \code{'observed'}.
#'@param standardise Logical, \code{FALSE} by default. If \code{TRUE}, report
#'the effect of a one standard deviation change in process c as a change in r's
#'own standard deviations, taking both scales from the stationary within
#'subject covariance \code{asymDIFFUSIONcov}. That makes effects between
#'processes on different scales comparable, and is what \code{\link{ctNetwork}}
#'reports by default. \code{FALSE} leaves the effects in the processes' own
#'units, so \code{impulseType='unit', standardise=FALSE} is the bare
#'\code{expm(DRIFT*t)}. A system with no stationary variance has nothing to
#'standardise by, and those draws are returned as NaN with a message. The plot
#'title states which scale was used either way.
#'@param cov Logical. If TRUE, covariances are returned instead of regression coefficients.
#'@param plot Logical. If TRUE, ggplots output using \code{\link{ctDiscreteParsPlot}}
#'instead of returning output.
#'@param cores Number of cpu cores to use for computing subject matrices.
#'If subject matrices were saved during fiting, not used.
#'@param state Evaluation point for any matrix cell that depends on the latent
#'state or a time dependent predictor -- see \code{\link{ctContextDependence}}.
#'\code{NULL} (the default) or \code{'T0MEANS'} uses the population initial
#'state, \code{'mean'} the mean smoothed latent state, \code{'asymptotic'} the
#'system's own fixed point; or supply a numeric state vector. Requires
#'\code{backend='julia'} and \code{subjects='popmean'}. For a linear model the
#'argument makes no difference, since no cell depends on the state.
#'@param method How the interval regressions are obtained. \code{'linearise'}
#'(the default) exponentiates the DRIFT matrix, \code{expm(DRIFT*t)} -- exact
#'for a linear model, and for a model whose DRIFT depends on the state, the
#'response of the model linearised at \code{state}. \code{'simulate'}
#'integrates the nonlinear system from \code{state} and from \code{state} plus
#'one unit on each process and differences the trajectories, which is the actual
#'model implied regression; it starts at 1 and decays to 0 for a stable process,
#'and reduces exactly to \code{'linearise'} when no cell depends on the state.
#'\code{'simulate'} requires \code{backend='julia'} and is far slower, so it
#'uses \code{nsamples} posterior draws with a small default.
#'@param ... additional plotting arguments to control \code{\link{ctDiscreteParsPlot}}
#'@examples
#' data.table::setDTthreads(1) #ignore this line
#' ctDiscretePars(ctstantestfit,times=seq(.5,4,.1),
#'  plot=TRUE,indices='CR')
#'
#'#modify plot
#'require(ggplot2)
#'g=ctDiscretePars(ctstantestfit,times=seq(.5,4,.1),
#'  plot=TRUE,indices='CR')
#'g= g+ labs(title='Cross effects')
#'print(g)
#'@details
#'Two of these combinations are quantities you may already know under other
#'names. \code{impulseType='unit'} with \code{standardise=FALSE} is the
#'discrete time autoregression and cross-lagged matrix, \code{expm(DRIFT*t)}.
#'\code{impulseType='observed'} with \code{standardise=TRUE} is the model
#'implied \emph{cross-correlation function},
#'\code{Cor(x_r(t+u), x_c(t))} -- the model's counterpart to the empirical
#'autocorrelations \code{\link{ctACF}} computes from the data, and equal to
#'the latent correlation matrix at \code{t=0}.
#'
#'If plot=TRUE, the function will return a ggplot2 object
#'(and hence needs to be printed if intended to display within a loop).
#'This can be modified by the various ggplot2 functions, or displayed using print(x).
#'@return With \code{plot = FALSE}, an array of the discrete time
#'coefficients, by sample, subject, time interval, row and column, of class
#'\code{ctDiscretePars}: \code{summary()} gives each curve's peak, the lag of
#'the peak, the lag at which it has fallen to half the peak (a half-life), and
#'the area under it, with their uncertainty (\code{\link{summary.ctDiscretePars}}).
#'With \code{plot = TRUE}, a ggplot object.
#'@aliases ctStanDiscretePars
#'@export
ctDiscretePars<-function(fit, subjects='popmean',
  times=seq(from=0,to=10,by=.1),
  nsamples=200,impulseType='unit',standardise=FALSE,
  cov=FALSE, plot=FALSE,cores=getOption("mc.cores", 2L),state=NULL,method='linearise',realid=TRUE,...,
  ctstanfitobj, observational){

  if(missing(fit)){
    if(missing(ctstanfitobj)) stop('fit must be supplied')
    warning('ctstanfitobj argument is deprecated, use fit instead')
    fit <- ctstanfitobj
  } else if(!missing(ctstanfitobj)) {
    stop('Use only one of fit or deprecated ctstanfitobj')
  }

  # `observational` said what the *design* was; `impulseType` says what the
  # companion matrix does, which is what actually varies and what the two
  # non-logical values never fitted into. The logical still resolves.
  if(!missing(observational)){
    if(!missing(impulseType)) stop(call.=FALSE,
      'Use only one of impulseType or deprecated observational')
    warning(call.=FALSE, "observational is deprecated, use impulseType='",
      .ctCompanionType(observational), "'")
    impulseType <- observational
  }
  impulseType <- .ctCompanionType(impulseType)

  method <- match.arg(method, c('linearise','simulate'))

  # method='simulate' answers a different question and takes a different route
  # to it: integrate the nonlinear system from a state and from that state plus
  # one unit, and difference. It reduces exactly to the linearised answer when
  # no cell depends on the state, so it is safe to ask for either way.
  if(method %in% 'simulate'){
    if(!'popmean' %in% subjects) stop(call.=FALSE,
      "method='simulate' applies to subjects='popmean' only.")
    ctmS <- .ctFitModelObject(fit)
    if(!ctmS$continuoustime) stop(call.=FALSE,
      "method='simulate' is for continuous time models.")
    if(standardise) .ctContextStandardiseMessage(fit)
    out <- .ctDiscreteParsSimulate(fit, times=times,
      state=if(is.null(state)) 'asymptotic' else state, nsamples=nsamples,
      impulseType=impulseType, standardise=standardise)
    times <- attr(out,'times')
    dimnames(out) <- list(Sample=seq_len(dim(out)[1]), Subject='popmean',
      `Time interval`=times, row=ctmS$latentNames, col=ctmS$latentNames)
    attributes(out)$impulseType <- impulseType
    attributes(out)$cov <- FALSE
    attributes(out)$standardise <- standardise
    attributes(out)$method <- 'simulate'
    out <- .ctContextAttach(out, fit)
    class(out) <- c('ctDiscretePars', 'array')
    if(plot) out <- ctDiscreteParsPlot(out, ...)
    return(out)
  }

  # `fit` may be a ctStanFit or a ctJuliaFit. Everything this
  # function needs is either the pop_* arrays that ctExtract() now returns for
  # all three, or model metadata; .ctFitModelObject() supplies the latter so the
  # body no longer reaches into $ctstanmodel and $standata directly.
  ctm <- .ctFitModelObject(fit)
  if(!ctm$continuoustime) times <- unique(round(times))
  type='discreteDRIFT'
  collapseSubjects=TRUE #consider this for a switch

  subjectLabels <- subjects
  if(!'popmean' %in% subjects){
    subjects <- .ctResolveSubjects(fit, subjects, realid)
    # ctExtract returns subjects in ctsem's order, so the labels follow it.
    ordering <- order(subjects)
    subjectLabels <- attr(subjects, 'ids')[ordering]
    subjects <- as.vector(subjects)[ordering]
  }

  extractSubjects <- subjects
  if('popmean' %in% extractSubjects) extractSubjects <- 'all'
  # An explicit evaluation point only means something for the population
  # matrices: a subject's are the ones its own filter pass ended with.
  stateArgs <- list()
  if(!is.null(state)){
    if(!inherits(fit,'ctJuliaFit')) stop(call.=FALSE,
      paste0("state= requires backend='julia'. ",
        "The stan model materialises its matrices during the fit, at T0MEANS, ",
        "and cannot be asked for another point afterwards."))
    if(!'popmean' %in% subjects) stop(call.=FALSE,
      paste0("state= applies to subjects='popmean' only; a subject's matrices ",
        "come from its own filter pass."))
    stateArgs$state <- .ctResolveState(fit,state)$state
  }
  e<-do.call(ctExtract,c(list(fit,subjectMatrices = subjects[1]!='popmean',cores=cores,
    nsamples = min(nsamples,.ctFitNsamples(fit)),
    subjects=extractSubjects, realid=FALSE),stateArgs))

  # A fit with no stationary variance to standardise by says so and carries on,
  # rather than returning a block of NaN. Same fallback as ctNetwork(), same
  # wording.
  if(standardise && is.null(e$pop_asymDIFFUSIONcov)){
    message('standardise=TRUE needs a stationary variance to standardise by, ',
      'and this fit has none; reporting unstandardised regressions.')
    standardise <- FALSE
  }

  nsubjects <- dim(e$subj_DRIFT)[2]
  if(is.null(nsubjects)) nsubjects=1
  if('popmean' %in% subjects) subjects='popmean'

  if(is.null(e$subj_DRIFT) && any(!subjects %in% 'popmean')) stop('No individual variation in DRIFT matrix found?? Try subjects="popmean"')

  niter=dim(e$pop_DRIFT)[1]
  latentNames=ctm$latentNames
  nlatent=length(latentNames)

  if(nsamples > niter) nsamples <- niter

  out<-list()

  #get all ctparameter matrices at once and remove unneeded subjects
  ctpars <- list()

  samples <- sample(1:dim(e[['pop_DRIFT']])[1],nsamples)

  # browser()
  for(matname in c('DRIFT','DIFFUSIONcov','asymDIFFUSIONcov')){ #,'CINT','T0MEANS', 'T0VAR','MANIFESTMEANS',if(!is.null(e$MANIFESTVAR)) 'MANIFESTVAR','LAMBDA', if(!is.null(e$TDPREDEFFECT)) 'TDPREDEFFECT')){
    if('popmean' %in% subjects || is.null(e[[paste0('subj_',matname)]])){
      ctpars[[matname]] <- e[[paste0('pop_',matname)]][samples,,,drop=FALSE]
    } else {
      if(dim(e[[paste0('subj_',matname)]])[2] != .ctFitNsubjects(fit)){ #if we computed subject parameters for only the specified subjects
        parsubjects <- 1:length(subjects)
      } else parsubjects <- subjects
      ctpars[[matname]] <- e[[paste0('subj_',matname)]][samples,parsubjects,,,drop=FALSE]
      # ctpars[[matname]]<-  array(ctpars[[matname]],dim=c(prod(dim(ctpars[[matname]])[1:2]),dim(ctpars[[matname]])[-1:-2]))
    }
    # ctpars[[matname]] <- ctpars[[matname]][sample(1:dim(ctpars[[matname]])[1],nsamples),,,drop=FALSE]
  }


  # Which evaluation point the DRIFT being exponentiated below actually came
  # from -- the population pass, or each subject's own filter pass. They are
  # different points and the reader is told which one they got.
  #
  # Kept in a variable and attached as well as messaged, because
  # ctDiscreteParsPlot() is a separate call: the producer was careful and the
  # figure it produced still left the session with no record of the point.
  stateLabel <- if(!'popmean' %in% subjects) .ctContextSubjectLabel else
    if(is.null(state)) .ctContextPopLabel else .ctResolveState(fit,state)$label
  .ctContextMessage(fit, stateLabel,
    paste0("expm(DRIFT*t) is therefore the transition of the model linearised ",
      "there, not the nonlinear system's own interval regression."))
  if(standardise) .ctContextStandardiseMessage(fit, stateLabel)

  out <- ctDiscreteParsDrift(ctpars,times, impulseType, standardise, cov=cov,discreteInput = ctm$continuoustime==FALSE)

  dimnames(out)<- list(Sample=samples, Subject=if('popmean' %in% subjects) subjects else subjectLabels,
    `Time interval`=times, row=latentNames, col=latentNames)

  attributes(out)$impulseType <- impulseType
  attributes(out)$cov <- cov
  attributes(out)$standardise <- standardise
  attributes(out)$method <- 'linearise'
  class(out) <- c('ctDiscretePars', 'array')
  attributes(out)$stateLabel <- stateLabel
  out <- .ctContextAttach(out, fit)

  if(plot) {

    out <- ctDiscreteParsPlot(out,
      # latentNames=ctstanfitobj$ctstanmodel$latentNames,
      ...)
  }
  return(out)
}

#' @export
ctStanDiscretePars <- ctDiscretePars



#' Summarise discrete time effect curves
#'
#' For each curve \code{\link{ctDiscretePars}} returned -- one per subject and
#' pair of processes -- the quantities read off it, per draw, and their
#' distribution over draws.
#'
#' @param object Output of \code{\link{ctDiscretePars}} with \code{plot = FALSE}.
#' @param probs Quantiles of each quantity over draws.
#' @param ... Not used.
#' @return A data frame, one row per subject, pair of processes (\code{row} the
#'   one affected, \code{col} the one acting) and quantity, with the mean, sd
#'   and quantiles over draws. The quantities:
#'   \describe{
#'     \item{\code{peak}}{The value of largest magnitude, with its sign.}
#'     \item{\code{peaktime}}{The lag at which it occurs.}
#'     \item{\code{halftime}}{The first lag after the peak at which the curve's
#'       magnitude has fallen to half the peak's, interpolated between lags.
#'       For an autoregression, which peaks at lag zero, its half-life: with
#'       \code{impulseType = 'observed', standardise = TRUE}, that of the
#'       autocorrelation. \code{NA} for a draw whose curve does not fall that
#'       far within the lags computed; \code{beyond} is the share of such
#'       draws, and the summary is over the others.}
#'     \item{\code{area}}{The area under the curve over the lags computed: for
#'       a unit impulse over long enough lags, the cumulative effect.}
#'   }
#'   The lags are those \code{ctDiscretePars} was given, so their spacing
#'   bounds the precision of \code{peaktime} and \code{halftime}, and their
#'   range what \code{halftime} and \code{area} can see.
#' @examples
#' dp <- ctDiscretePars(ctstantestfit, times = seq(0, 5, .1), nsamples = 20)
#' summary(dp)
#' @export
summary.ctDiscretePars <- function(object, probs = c(.025, .5, .975), ...) {
  times <- as.numeric(dimnames(object)[[3L]])
  dims <- dim(object)
  names <- dimnames(object)
  rows <- list()
  for (i in seq_len(dims[2L])) for (r in seq_len(dims[4L])) for (c in seq_len(dims[5L])) {
    V <- matrix(unclass(object)[, i, , r, c], nrow = dims[1L])
    quantities <- .ctCurveQuantities(V, times)
    for (q in names(quantities)) {
      v <- quantities[[q]]
      ok <- is.finite(v)
      summarised <- if (any(ok)) c(mean(v[ok]), if (sum(ok) > 1L) stats::sd(v[ok]) else
        NA_real_, stats::quantile(v[ok], probs, names = FALSE)) else
        rep(NA_real_, 2L + length(probs))
      rows[[length(rows) + 1L]] <- data.frame(Subject = names[[2L]][i],
        row = names[[4L]][r], col = names[[5L]][c], quantity = q,
        t(summarised), beyond = if (q == "halftime") mean(!ok) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  colnames(out)[5L:(6L + length(probs))] <- c("mean", "sd",
    paste0(100 * probs, "%"))
  rownames(out) <- NULL
  out
}

# Each draw's peak, its lag, the lag the curve falls to half the peak after
# it, and the area under it: one curve per row of `V`, over `times`.
.ctCurveQuantities <- function(V, times) {
  A <- abs(V)
  peakcol <- max.col(A, ties.method = "first")
  index <- cbind(seq_len(nrow(V)), peakcol)
  peak <- V[index]
  half <- abs(peak) / 2
  after <- A <= half & col(A) > peakcol
  first <- max.col(after, ties.method = "first")
  reached <- rowSums(after) > 0 & abs(peak) > 0
  halftime <- rep(NA_real_, nrow(V))
  for (s in which(reached)) {
    j <- first[s]
    a0 <- A[s, j - 1L]
    a1 <- A[s, j]
    halftime[s] <- times[j - 1L] + (a0 - half[s]) / (a0 - a1) *
      (times[j] - times[j - 1L])
  }
  steps <- diff(times)
  area <- if (length(times) > 1L) as.numeric(((V[, -1L, drop = FALSE] +
    V[, -ncol(V), drop = FALSE]) / 2) %*% steps) else rep(0, nrow(V))
  list(peak = peak, peaktime = times[peakcol], halftime = halftime, area = area)
}

#' @export
print.ctDiscretePars <- function(x, ...) {
  y <- x
  class(y) <- NULL
  print(y, ...)
  invisible(x)
}

# Matrices batched over their first dimension: an array [n, row, col] holds n
# matrices, and these act on all of them at once, so ctDiscreteParsDrift()
# loops over lags and the matrices' own few rows and columns, never over draws.
.ctBatchMul <- function(A, B) {
  nl <- dim(A)[2]
  out <- array(0, dim(A))
  for (i in seq_len(nl)) for (j in seq_len(nl)) {
    s <- A[, i, 1] * B[, 1, j]
    if (nl > 1) for (k in 2:nl) s <- s + A[, i, k] * B[, k, j]
    out[, i, j] <- s
  }
  out
}

.ctBatchIdentity <- function(n, nl) {
  out <- array(0, c(n, nl, nl))
  for (i in seq_len(nl)) out[, i, i] <- 1
  out
}

# expm() of every matrix: scaled by a power of two until every 1-norm is under
# one half, a Taylor series of order 16 (remainder under 1e-19 there), then
# squared back.
.ctBatchExpm <- function(A) {
  n <- dim(A)[1]
  nl <- dim(A)[2]
  norm <- numeric(n)
  for (j in seq_len(nl)) {
    s <- 0
    for (i in seq_len(nl)) s <- s + abs(A[, i, j])
    norm <- pmax(norm, s)
  }
  biggest <- max(c(0, norm[is.finite(norm)]))
  squarings <- if (biggest > 0.5) ceiling(log2(biggest / 0.5)) else 0
  A <- A / 2^squarings
  I <- .ctBatchIdentity(n, nl)
  X <- I
  for (k in 16:1) X <- I + .ctBatchMul(A, X) / k
  for (q in seq_len(squarings)) X <- .ctBatchMul(X, X)
  X
}

# The p-th power of every matrix, p a whole number, by repeated squaring.
.ctBatchPower <- function(A, p) {
  out <- .ctBatchIdentity(dim(A)[1], dim(A)[2])
  while (p > 0) {
    if (p %% 2 == 1) out <- .ctBatchMul(out, A)
    p <- p %/% 2
    if (p > 0) A <- .ctBatchMul(A, A)
  }
  out
}

ctDiscreteParsDrift<-function(ctpars,times, impulseType,  standardise,cov=FALSE,
  types='dtDRIFT',discreteInput=FALSE, quiet=FALSE){

  nl=dim(ctpars$DRIFT)[3]

  if(!quiet) message('Computing temporal regression coefficients for ', dim(ctpars$DRIFT)[1],' samples')

  lapply(names(ctpars),function(x){ #add in extra dim if only 3 dims (e.g. when not individually varying)
    dm=dim(ctpars[[x]])
    if(length(dm)==3){
      ctpars[[x]] <<- array(ctpars[[x]],dim=c(dm[1],1,dm[2:3]))
    }
  })

  nsubs <- lapply(ctpars,function(x) dim(x)[2])


  nonstationary <- 0L
  companionType <- .ctCompanionType(impulseType)

  if('dtDRIFT' %in% types){
    nt <- length(times)
    S <- dim(ctpars$DRIFT)[1]
    J <- max(unlist(nsubs))
    N <- S * J
    dt <- array(NA_real_, dim=c(S, J, nt, nl, nl))
    # Every draw of every subject at once, as [draw x subject, row, col], a
    # subject past a matrix's own taking its last, as before. A loop over
    # draw, subject and lag with an expm() each took 20 minutes for 10,000
    # draws of three subjects at 361 lags.
    flat <- function(m) array(m[, pmin(seq_len(J), dim(m)[2]), , , drop = FALSE], c(N, nl, nl))
    D <- flat(ctpars$DRIFT)
    bad <- logical(N)

    C <- NULL
    if(!identical(companionType,'unit')){
      # Every variant is dtDRIFT %*% C; only the companion matrix C
      # differs. See R/ctCompanionShock.R for what each one means.
      Q <- flat(ctpars$DIFFUSIONcov)
      A <- flat(ctpars$asymDIFFUSIONcov)
      if(companionType %in% c('observed', 'shock')){
        # .ctCovarianceRegression() of every draw at once: each column of the
        # covariance over its variance, 0 where that is not a positive number,
        # ones on the diagonal; 'observed' refuses a draw whose stationary
        # variances are not all finite and non-negative, as .ctCompanionMatrix()
        # does.
        M <- if(companionType == 'observed') A else Q
        variances <- matrix(NA_real_, N, nl)
        for(r in seq_len(nl)) variances[, r] <- M[, r, r]
        if(companionType == 'observed')
          bad <- bad | apply(!is.finite(variances) | variances < 0, 1, any)
        inverse <- 1 / variances
        inverse[!is.finite(variances) | variances <= 0] <- 0
        C <- M
        for(cc in seq_len(nl)) C[, , cc] <- M[, , cc] * inverse[, cc]
        C[!is.finite(C)] <- 0
        for(r in seq_len(nl)) C[, r, r] <- 1
      } else {
        C <- array(NaN, c(N, nl, nl))
        for(n in seq_len(N)){
          Cn <- .ctCompanionMatrix(companionType, matrix(Q[n,,],nl,nl), matrix(A[n,,],nl,nl), nl)
          if(is.null(Cn)) bad[n] <- TRUE else C[n,,] <- Cn
        }
      }
    }

    if(standardise) {
      # Response of r per one sd_c change in c, in sd_r units: S^-1 M S,
      # with the sds of the processes themselves (asymDIFFUSION). Applied
      # after the companion step so the two compose.
      #
      # A frozen DRIFT can be non-stationary at the point it was frozen at
      # while the nonlinear system it came from is perfectly well behaved
      # -- often the point of the nonlinearity -- so this reports NaN for
      # the affected draw and names the likely cause rather than aborting
      # and blaming the model, as julia's _ctsem_asymptotics already does.
      A <- flat(ctpars$asymDIFFUSIONcov)
      variances <- matrix(NA_real_, N, nl)
      for(r in seq_len(nl)) variances[, r] <- A[, r, r]
      bad <- bad | apply(variances < 0, 1, function(z) isTRUE(any(z)))
      sdv <- suppressWarnings(sqrt(variances + 1e-10))
      scale <- array(NA_real_, c(N, nl, nl))
      for(r in seq_len(nl)) for(cc in seq_len(nl)) scale[, r, cc] <- sdv[, cc] / sdv[, r]
    }
    nonstationary <- nonstationary + nt * sum(bad)

    # The lags in increasing order, each the last times expm(DRIFT * step),
    # or DRIFT to the step's power in discrete time. A step of the same length
    # as one before, to 13 figures, reuses it: an evenly spaced grid computes
    # one or two.
    E <- .ctBatchIdentity(N, nl)
    last <- 0
    steps <- list()
    for(ti in order(times)){
      step <- times[ti] - last
      key <- format(signif(step, 13), digits = 15)
      X <- steps[[key]]
      if(is.null(X)){
        X <- if(discreteInput) .ctBatchPower(D, round(step)) else .ctBatchExpm(D * step)
        if(length(steps) < 8) steps[[key]] <- X
      }
      E <- .ctBatchMul(E, X)
      last <- times[ti]
      Et <- E
      if(!is.null(C)) Et <- .ctBatchMul(Et, C)
      if(standardise) Et <- Et * scale
      if(cov) {
        Ec <- Et
        for(r in seq_len(nl)) for(q in seq_len(nl)) Ec[, r, q] <- rowSums(Et[, r, , drop=FALSE] * Et[, q, , drop=FALSE])
        Et <- Ec
      }
      Et[bad, , ] <- NaN
      dt[, , ti, , ] <- Et
    }
    ctpars$dtDRIFT <- dt
  } #end dtdrift

  if(nonstationary > 0 && !quiet) message(
    nonstationary, ' of ', length(ctpars$dtDRIFT)/prod(dim(ctpars$DRIFT)[3:4]),
    ' standardisations returned NaN: the asymptotic diffusion has negative ',
    'diagonals, so there is no stationary variance to standardise by. For a ',
    'model whose DRIFT depends on the state this usually reflects the point it ',
    'was linearised at rather than the model -- try another state, or ',
    'standardise=FALSE.')

  return(ctpars$dtDRIFT)
}

#'ctDiscreteParsPlot
#'
#'Plots the output from \code{\link{ctDiscretePars}}, for model implied regression strengths at specified times for continuous time models fit with ctStanFit.
#'\code{ctStanDiscreteParsPlot} is maintained as a backward-compatible alias.
#'
#'@param x list object returned from \code{\link{ctDiscretePars}}.
#'@param indices Either a string specifying type of plot to create, or an n by 2
#'matrix specifying which indices of the output matrix to plot.
#''AR' specifies all diagonals, for discrete time autoregression parameters.
#''CR' specifies all off-diagonals,for discrete time cross regression parameters.
#''all' plots all AR and CR effects at once.
#'@param quantiles numeric vector of length 3, with values between 0 and 1, specifying which quantiles to plot.
#'The default plots 95\% credible intervals and the posterior median at 50\%.
#'@param latentNames Vector of character strings denoting names for the latent variables.
#''auto' just uses eta1 eta2 etc.
#'@param polygonalpha Numeric between 0 and 1 to multiply the alpha of
#'the fill.
#'@param ylab y label.
#'@param xlab x label.
#'@param ylim Two values, the lower and upper limits of the y axis. The plot
#'zooms to them; intervals reaching past them are drawn to the edge.
#'@param facets May be 'Subject' or 'Effect'.
#'@param colour Character string denoting how colour varies. 'Effect' or 'Subject'.
#'@param title Character string. 'auto' generates automatically, NULL can be used to disable title.
#'@param splitSubjects if TRUE, subjects are plotted separately, if FALSE they are combined.
#'@param ggcode if TRUE, returns a list containing the data.table to plot (for each curve and time interval,
#'the three \code{quantiles} of the draws, as \code{ymin}, \code{value} and \code{ymax}), and a character string that can be
#'evaluated (with the necessary arguments such as ylab etc filled in). For modifying plots.
#'@return A ggplot2 object. This can be modified by the various ggplot2 functions, or displayed using print(x).
#'@examples
#' data.table::setDTthreads(1) #ignore this line
#'x <- ctDiscretePars(ctstantestfit)
#'ctDiscreteParsPlot(x, indices='CR')
#'
#'#to modify plot:
#'g <- ctDiscreteParsPlot(x, indices='CR') +
#'  ggplot2::labs(title='My ggplot modification')
#'print(g)
#'
#'@aliases ctStanDiscreteParsPlot
#'@export

ctDiscreteParsPlot<- function(x,indices='all',
  quantiles=c(.025,.5,.975), latentNames='auto',
  ylab='Coefficient',xlab='Time interval',ylim=NA,facets=NA,splitSubjects=TRUE,
  colour='Effect',title='auto',
  polygonalpha=.1,ggcode=FALSE){

  if(is.data.frame(indices)) indices <- as.matrix(indices)

  impulseLabel <- NULL
  if(!is.null(title)){
    if(title %in% 'auto'){
      # Whether the curves are in standard deviations or in the processes' own
      # units is not recoverable from the picture, and the two differ by a
      # factor per process pair. NULL for an object made before the attribute
      # existed, where the scale is unknown and is better not claimed.
      # `observational` is read as a fallback for the same reason.
      std <- attributes(x)$standardise
      impulse <- attributes(x)$impulseType
      if(is.null(impulse)) impulse <- attributes(x)$observational
      # An object old enough to carry neither is the default, which is what it
      # would have been computed with.
      if(is.null(impulse)) impulse <- 'unit'
      # The impulse goes in the subtitle: in the title with the scale, the
      # longest ('observed change of 1 in one process, others follow') was
      # cut off at the 480-pixel width of a default png().
      title= paste0('Temporal ',ifelse(attributes(x)$cov,'covariance','regressions'),
        if(isTRUE(std)) ', standardised' else
          if(isFALSE(std)) ', unstandardised' else '')
      impulseLabel <- .ctCompanionLabel(impulse, std)
    }
  }else title=title

  # A state dependent object's figure has to carry its evaluation point. These
  # curves are expm() of a DRIFT that was linearised somewhere, and which
  # somewhere is not recoverable from the picture; ctDiscretePars() attached the
  # cells and the label for exactly this and nothing downstream read them.
  # Silent for a linear model, whose DRIFT is the same everywhere.
  .cells <- attributes(x)$contextDependent
  .conditional <- !is.null(.cells) && nrow(.cells) > 0 &&
    any(.cells$kind %in% .ctContextProblemKinds)
  subtitle <- if(isTRUE(.conditional) && !is.null(attributes(x)$stateLabel))
    paste0('State dependent ',
      paste0(.ctContextReportableMatrices(
        .cells[.cells$kind %in% .ctContextProblemKinds, , drop = FALSE]),
        collapse = '/'),
      ': linearised at ', attributes(x)$stateLabel, '.') else NULL

  nlatent=dim(x)[5]

  if(latentNames[1]=='auto') latentNames=dimnames(x)$row
  dimnames(x)$row <- dimnames(x)$col <- latentNames

  if(all(indices=='AR')) indices <- matrix(1:nlatent,nrow=nlatent,ncol=2)

  if(all(indices=='CR')) indices <- cbind(
    rep(1:nlatent,nlatent)[-seq(1,nlatent^2,nlatent+1)],
    rep(1:nlatent,each=nlatent-1))

  if(indices[1]=='all') indices <- cbind(
    rep(1:nlatent,nlatent),
    rep(1:nlatent,each=nlatent))

  if(!'array' %in% class(indices) && !'matrix' %in% class(indices)){#interpret as individual columns
    indices <- cbind(
      rep(1:nlatent,length(unique(indices))),
      rep(unique(indices),each=nlatent))
  }


  # One row per curve and time interval, holding the quantiles of the draws
  # there: what the plot draws. The draws themselves were once its data, and
  # every print summarised them again -- 10.8 million rows for 10,000 draws of
  # three subjects at 361 intervals. Non-finite draws are left out, as ggplot
  # left them out of the summary it made.
  levs <- levels(factor(latentNames))
  times <- as.numeric(dimnames(x)$`Time interval`)
  subjects <- if(splitSubjects) as.character(dimnames(x)$Subject) else '1'
  pairs <- unique(matrix(indices[, 1:2], ncol = 2))
  ym <- rbindlist(lapply(seq_len(nrow(pairs)), function(k) {
    v <- x[, , , pairs[k, 1], pairs[k, 2], drop = FALSE]
    dim(v) <- dim(v)[1:3] # sample, subject, time
    if(!splitSubjects) dim(v) <- c(dim(v)[1] * dim(v)[2], 1L, dim(v)[3])
    q <- apply(v, c(2, 3), function(z) {
      z <- z[is.finite(z)]
      if(length(z)) stats::quantile(z, quantiles[1:3], names = FALSE) else rep(NA_real_, 3)
    })
    data.table(Subject = rep(subjects, times = length(times)),
      `Time interval` = rep(times, each = length(subjects)),
      row = latentNames[pairs[k, 1]], col = latentNames[pairs[k, 2]],
      ymin = c(q[1, , ]), value = c(q[2, , ]), ymax = c(q[3, , ]))
  }))
  ym$row <- factor(ym$row, levels = levs)
  ym$col <- factor(ym$col, levels = levs)
  ym$Effect <- interaction(ym$row,ym$col)
  ym$Subject <- factor(ym$Subject)

  # The interval first and the median line second: code that adjusts the
  # layers by position, as the ctsem book does, relies on that order.
  g<-paste0('ggplot2::ggplot(data = ym,mapping=aes(x=`Time interval`,
    colour=Effect,group=interaction(Effect, Subject),
    fill=Effect))+
    theme_bw()+
    ylab("',ylab,'")+
    xlab("',xlab,'")+
    ggplot2::labs(title = "',title,'")+
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
      alpha= ',polygonalpha,',
      linetype=3)+
    geom_line(aes(y = value))')

  #needs data input as well...
  # if(rug) g <- paste0(g,'+ geom_rug(data= as.data.table(hist(data.table(ym)[order(Subject,`Time interval`),][,.(timeDiff=c(diff(`Time interval`))),by=Subject]$timeDiff,plot=F)[2:4])[counts > 0,],
  #  aes(x = mids,alpha=density,size=density/max(density)),inherit.aes=F,linetype=1,show.legend = FALSE)')

  # Each a vector: a test of `!is.na()` on one stopped every call that gave
  # two limits. The limits zoom (`coord_cartesian`) rather than set the scale
  # (`ylim()`), which drops every interval that leaves them and so cut the band.
  if(!all(is.na(facets))) g <- paste0(g,'+ facet_wrap(facets)')

  if(!all(is.na(ylim))) g <- paste0(g,' + ggplot2::coord_cartesian(ylim = ylim)')

  # The impulse, and under it any evaluation point.
  subtitle <- c(impulseLabel, subtitle)
  if(length(subtitle) && !is.null(title)) g <- paste0(g,
    '+ ggplot2::labs(subtitle = "', paste(subtitle, collapse = '\\n'), '")')

  if(!ggcode) g <- eval(parse(text=g)) else g <- list(dt=ym,ggcode=g)

  return(g)
}

#' @export
ctStanDiscreteParsPlot <- ctDiscreteParsPlot
# } else {

#
#   #blank plot
#   rm(plot) #otherwise can't find function
#   blankargs=plotcontrol
#   blankargs$xlim=range(times)
#   blankargs$y=NA
#   blankargs$x=NA
#   do.call(plot,blankargs)
#   if(grid) {
#     grid()
#     par(new=TRUE)
#     do.call(plot,blankargs)
#     par(new=FALSE)
#   }
#
#   ####plotting confidence region
#   if(polygon) {
#     cc=0
#     ccup=TRUE
#     for(indexi in c(1:nrow(indices))){
#       cc=ifelse(ccup,cc+1,cc-1)
#       if(indexi==nrow(indices)) ccup=FALSE
#       ri=indices[indexi,1]
#       ci=indices[indexi,2]
#       polygonargs<-polygoncontrol
#       polygonargs$x=times
#       polygonargs$y=x[ri,ci,,2]
#       polygonargs$ylow=x[ri,ci,,1]
#       polygonargs$yhigh=x[ri,ci,,length(quantiles)]
#       polygonargs$col=grDevices::adjustcolor(colvec[cc],alpha.f=max(c(.004,polygonalpha/(2*sqrt(polygonargs$steps)))))
#       do.call(ctPoly,polygonargs)
#     }
#   }
#
#   ####plotting quantile lines
#   for(qi in 1:3){
#     cc=0
#     for(indexi in 1:nrow(indices)){
#       cc=cc+1
#       ri=indices[indexi,1]
#       ci=indices[indexi,2]
#
#       plotargs<-plotcontrol
#       plotargs$x=times
#       plotargs$y=x[ri,ci,,qi]
#       plotargs$lty=ltyvec[cc]
#       plotargs$col=ifelse(qi!= 2,grDevices::adjustcolor(colvec[cc],alpha.f=.5) ,colvec[cc])
#       plotargs$lwd=ifelse(qi!= 2,1, lwdvec[cc])
#       do.call(points,plotargs)
#     }}
#
#
#
#   legendcontrol$legend=paste0(latentNames[indices[,1]],'_',latentNames[indices[,2]])
#   legendcontrol$text.col=colvec
#   legendcontrol$col=colvec
#   legendcontrol$lty = ltyvec
#   legendcontrol$lwd=lwdvec
#
#   if(legend) do.call(graphics::legend,legendcontrol)
# }#end if not gg
# }
