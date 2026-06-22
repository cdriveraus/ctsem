ctEBnum <- function(x){
  sprintf('%.17g', x)
}

ctEBrawParnames <- function(fit){
  pnames <- try(getparnames(fit), silent=TRUE)
  if('try-error' %in% class(pnames) || length(pnames) != length(fit$stanfit$rawest)){
    pnames <- paste0('par', seq_along(fit$stanfit$rawest))
  }
  pnames
}

ctEBfreeRows <- function(pars, parnames){
  which(
    !is.na(pars$param) &
      is.na(pars$value) &
      pars$param %in% parnames &
      !grepl('[', pars$param, fixed=TRUE))
}

ctEBbaseTransform <- function(transform){
  if(length(transform) != 1 || is.na(transform)) return(NA_character_)
  transform <- as.character(transform)
  transformint <- suppressWarnings(as.integer(transform))
  if(!is.na(transformint)) {
    transform <- tformshapes(singletext=TRUE, transform=transformint)
  }
  transform
}

ctEBshiftTransform <- function(transform, rawmean, rawsd){
  transform <- ctEBbaseTransform(transform)
  if(is.na(transform)) return(NA_character_)
  rawexpr <- paste0('(', ctEBnum(rawmean), '+', ctEBnum(rawsd), '*param)')
  out <- gsub('\\bparam\\b', rawexpr, transform, perl=TRUE)
  if(identical(out, transform)) {
    warning('Could not find param in transform: ', transform)
  }
  out <- try(Simplify(out), silent=TRUE)
  if('try-error' %in% class(out)) out <- gsub('\\bparam\\b', rawexpr, transform, perl=TRUE)
  as.character(out)
}

ctEBadjustModel <- function(model, rawstats, sdscale=c('unit','rawsd'), minsd=1e-6){
  sdscale <- match.arg(sdscale)
  adjusted <- model
  rawstats <- rawstats[match(unique(rawstats$param), rawstats$param),,drop=FALSE]
  rownames(rawstats) <- rawstats$param
  rows <- ctEBfreeRows(adjusted$pars, rawstats$param)
  
  for(ri in rows){
    p <- adjusted$pars$param[ri]
    rmean <- rawstats[p, 'mean']
    rsd <- max(rawstats[p, 'sd'], minsd, na.rm=TRUE)
    adjusted$pars$transform[ri] <- ctEBshiftTransform(adjusted$pars$transform[ri], rmean, rsd)
    adjusted$pars$sdscale[ri] <- switch(sdscale,
      unit=1,
      rawsd=rsd)
  }
  
  adjusted$empiricalbayes <- list(
    rawstats=rawstats,
    sdscale=sdscale,
    minsd=minsd,
    note='Transforms evaluate the original parameter transform at rawmean + rawsd * param. EB covariance is reported but not represented in the adjusted model.')
  adjusted
}

ctEBrawMatrix <- function(fits, parnames, use=c('rawest','rawposterior')){
  use <- match.arg(use)
  
  if(use == 'rawest'){
    raw <- do.call(rbind, lapply(fits, function(fit) fit$stanfit$rawest))
    if(ncol(raw) != length(parnames)) stop('Raw point estimates do not match parnames')
    colnames(raw) <- parnames
    rownames(raw) <- names(fits)
    return(raw)
  }
  
  raw <- do.call(rbind, lapply(seq_along(fits), function(i){
    fit <- fits[[i]]
    samples <- ctStanRawSamples(fit)
    if(ncol(samples) != length(parnames)) {
      stop('Raw posterior samples do not match parnames for subject ', names(fits)[i])
    }
    colnames(samples) <- parnames
    rownames(samples) <- paste0(names(fits)[i], '.', seq_len(nrow(samples)))
    samples
  }))
  raw
}

ctEBrawStats <- function(raw, probs, location=c('mean','median'),
  scale=c('sd','mad','iqr')){
  location <- match.arg(location)
  scale <- match.arg(scale)
  centerfun <- switch(location,
    mean=function(x) mean(x, na.rm=TRUE),
    median=function(x) stats::median(x, na.rm=TRUE))
  scalefun <- switch(scale,
    sd=function(x) stats::sd(x, na.rm=TRUE),
    mad=function(x) stats::mad(x, center=stats::median(x, na.rm=TRUE),
      constant=1.4826, na.rm=TRUE),
    iqr=function(x) stats::IQR(x, na.rm=TRUE) / 1.349)
  out <- data.frame(
    param=colnames(raw),
    mean=as.numeric(apply(raw, 2, centerfun)),
    sd=as.numeric(apply(raw, 2, scalefun)),
    stringsAsFactors=FALSE)
  qs <- t(apply(raw, 2, stats::quantile, probs=probs, na.rm=TRUE))
  colnames(qs) <- paste0(probs * 100, '%')
  data.frame(out, qs, check.names=FALSE)
}

ctEBrobustRaw <- function(raw, outlierMAD=6, outlierQuantiles=c(.025,.975),
  winsorize=TRUE){
  cleaned <- raw
  report <- data.frame(
    param=colnames(raw),
    lower=NA_real_,
    upper=NA_real_,
    nlower=0L,
    nupper=0L,
    nchanged=0L,
    action=ifelse(winsorize, 'winsorized', 'removed'),
    stringsAsFactors=FALSE)
  
  for(j in seq_len(ncol(raw))){
    x <- raw[,j]
    finite <- is.finite(x)
    if(!any(finite)) next
    
    lower <- -Inf
    upper <- Inf
    
    if(!is.null(outlierQuantiles) && length(outlierQuantiles) == 2 &&
        all(is.finite(outlierQuantiles))){
      qbounds <- stats::quantile(x[finite], probs=outlierQuantiles,
        na.rm=TRUE, names=FALSE)
      lower <- max(lower, qbounds[1])
      upper <- min(upper, qbounds[2])
    }
    
    if(!is.null(outlierMAD) && is.finite(outlierMAD) && outlierMAD > 0){
      center <- stats::median(x[finite], na.rm=TRUE)
      madsd <- stats::mad(x[finite], center=center, constant=1.4826,
        na.rm=TRUE)
      if(!is.finite(madsd) || madsd <= 0) {
        madsd <- stats::IQR(x[finite], na.rm=TRUE) / 1.349
      }
      if(is.finite(madsd) && madsd > 0){
        lower <- max(lower, center - outlierMAD * madsd)
        upper <- min(upper, center + outlierMAD * madsd)
      }
    }
    
    low <- finite & x < lower
    high <- finite & x > upper
    report$lower[j] <- lower
    report$upper[j] <- upper
    report$nlower[j] <- sum(low)
    report$nupper[j] <- sum(high)
    report$nchanged[j] <- sum(low) + sum(high)
    
    if(report$nchanged[j] > 0){
      if(winsorize){
        cleaned[low,j] <- lower
        cleaned[high,j] <- upper
      } else {
        cleaned[low | high,j] <- NA_real_
      }
    }
  }
  
  list(raw=cleaned, report=report)
}

ctEBidentityRawMap <- function(parnames){
  data.frame(param=parnames, mean=0, sd=1, stringsAsFactors=FALSE)
}

ctEBmapRaw <- function(raw, rawmap){
  if(is.null(rawmap)) return(raw)
  mapped <- raw
  rawmap <- rawmap[match(unique(rawmap$param), rawmap$param),,drop=FALSE]
  rownames(rawmap) <- rawmap$param
  for(p in intersect(colnames(mapped), rawmap$param)){
    mapped[,p] <- rawmap[p, 'mean'] + rawmap[p, 'sd'] * mapped[,p]
  }
  mapped
}

ctEBpriorRawStats <- function(raw, ebRobust=TRUE, ebOutlierMAD=6,
  ebOutlierQuantiles=c(.025,.975), ebWinsorize=TRUE, minsd=1e-6,
  probs=c(.025,.5,.975)){
  
  rawForEB <- raw
  outliers <- NULL
  if(ebRobust){
    robustraw <- ctEBrobustRaw(raw, outlierMAD=ebOutlierMAD,
      outlierQuantiles=ebOutlierQuantiles, winsorize=ebWinsorize)
    rawForEB <- robustraw$raw
    outliers <- robustraw$report
  }
  rawstats <- ctEBrawStats(rawForEB, probs=probs,
    location=ifelse(ebRobust, 'median', 'mean'),
    scale=ifelse(ebRobust, 'mad', 'sd'))
  fallbackstats <- ctEBrawStats(rawForEB, probs=probs,
    location='median', scale='iqr')
  badsd <- !is.finite(rawstats$sd) | rawstats$sd < minsd
  rawstats$sd[badsd] <- fallbackstats$sd[badsd]
  rawstats$sd[!is.finite(rawstats$sd) | rawstats$sd < minsd] <- minsd
  
  list(raw=raw, rawForEB=rawForEB, rawstats=rawstats, outliers=outliers)
}

ctEBprogressReporter <- function(stage, total, enabled=TRUE){
  enabled <- isTRUE(enabled) && total > 0
  lastwidth <- 0L
  force(stage)
  force(total)
  function(done, finished=FALSE){
    if(!enabled) return(invisible(NULL))
    pct <- floor(100 * done / total)
    line <- sprintf('%s: %d/%d subjects (%d%%)', stage, done, total, pct)
    padding <- strrep(' ', max(0L, lastwidth - nchar(line)))
    cat('\r', line, padding, sep='', file=stderr())
    if(finished) cat('\n', file=stderr())
    lastwidth <<- nchar(line)
    invisible(NULL)
  }
}

ctEBfitSubjects <- function(subjects, datalong, subjectIDname, fitargs, cores=1,
  verbose=0, pass='', progress=TRUE){
  fitone <- function(subi){
    if(verbose > 1) {
      message('Fitting subject ', subi, if(nchar(pass)) paste0(' (', pass, ')') else '')
    }
    datasi <- datalong[datalong[[subjectIDname]] %in% subi,,drop=FALSE]
    do.call(ctFit, c(list(datalong=datasi), fitargs))
  }
  
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores)) cores <- 1L
  cores <- min(length(subjects), max(1L, cores))
  stage <- paste0('ctEmpiricalBayesFit', if(nchar(pass)) paste0(' ', pass) else '')
  progressfun <- ctEBprogressReporter(stage=stage, total=length(subjects),
    enabled=progress)
  progressdone <- FALSE
  done <- 0L
  on.exit({
    if(isTRUE(progress) && !progressdone) progressfun(done, finished=TRUE)
  }, add=TRUE)
  progressfun(0L)
  
  if(cores > 1){
    cl <- parallelly::makeClusterPSOCK(cores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(ctsem)))
    fits <- vector('list', length(subjects))
    nworkers <- min(cores, length(subjects))
    nextjob <- 1L
    for(worker in seq_len(nworkers)){
      parallel:::sendCall(cl[[worker]], fitone, list(subjects[[nextjob]]),
        tag=nextjob)
      nextjob <- nextjob + 1L
    }
    while(done < length(subjects)){
      result <- parallel:::recvOneResult(cl)
      done <- done + 1L
      fits[[result$tag]] <- result$value
      progressfun(done)
      if(nextjob <= length(subjects)){
        parallel:::sendCall(cl[[result$node]], fitone,
          list(subjects[[nextjob]]), tag=nextjob)
        nextjob <- nextjob + 1L
      }
    }
  } else {
    fits <- vector('list', length(subjects))
    for(i in seq_along(subjects)){
      fits[[i]] <- fitone(subjects[[i]])
      done <- i
      progressfun(done)
    }
  }
  progressfun(done, finished=TRUE)
  progressdone <- TRUE
  names(fits) <- as.character(subjects)
  fiterrors <- which(vapply(fits, inherits, logical(1), 'try-error'))
  if(length(fiterrors) > 0) {
    stop('Subject fit failed for subject ', names(fits)[fiterrors[1]], ': ',
      as.character(fits[[fiterrors[1]]]))
  }
  fits
}

ctEBfitArgsOptimDefaults <- function(fitargs, stochastic=FALSE,
  firstpass=FALSE){
  if(is.null(fitargs$optimcontrol)) fitargs$optimcontrol <- list()
  if(!is.list(fitargs$optimcontrol)) stop('optimcontrol must be a list')
  
  if(is.null(fitargs$optimcontrol$stochastic)) {
    fitargs$optimcontrol$stochastic <- stochastic
  }
  if(firstpass) fitargs$optimcontrol$estonly <- TRUE
  
  fitargs
}

#' Empirical Bayes subject-wise ctsem fits
#'
#' Fits one ctsem model per subject using the model prior, estimates the
#' empirical marginal distribution of the raw parameters, then fits each subject
#' again using the resulting empirical Bayes prior.
#'
#' @param datalong Long format data containing multiple subjects.
#' @param ctstanmodel Model object from \code{\link{ctModel}}. Time independent
#' predictors are not supported. A random-effect-free copy of this model is used
#' for the per-subject fits.
#' @param subjects Vector of subject identifiers to fit, or \code{'all'}.
#' @param priors Logical. Passed to \code{\link{ctFit}}; defaults to \code{TRUE}.
#' @param optimize Logical. Passed to \code{\link{ctFit}}; defaults to
#' \code{TRUE}.
#' @param cores Number of subjects to fit in parallel. Each individual
#' subject-level \code{\link{ctFit}} call uses one core.
#' @param subjectFitArgs Named list of additional arguments passed to each
#' \code{\link{ctFit}} call. For optimized fits, \code{optimcontrol$stochastic}
#' defaults to \code{FALSE} for all EB passes unless supplied here or in
#' \code{...}. First-pass fits force \code{optimcontrol$estonly=TRUE}.
#' @param Npasses Total number of subject-wise fitting passes. The default
#' \code{2} fits once with the model prior, builds a marginal empirical Bayes
#' prior, then fits once with that prior. Values above \code{2} repeatedly map
#' the previous pass estimates back to the original raw scale, rebuild the
#' marginal EB prior, and refit.
#' @param ebUse \code{'rawest'} to build the empirical Bayes prior from first
#' pass point estimates, or \code{'rawposterior'} to pool raw posterior samples.
#' @param ebRobust Logical. If TRUE, the empirical Bayes prior is built from
#' robust raw summaries after outlier handling.
#' @param ebOutlierMAD Positive numeric. Raw values farther than this
#' many MAD-scaled deviations from the median are treated as outliers. Use
#' \code{Inf} or \code{NULL} to disable this rule.
#' @param ebOutlierQuantiles Length two numeric vector of lower and upper
#' quantiles used to bound first-pass raw values, or \code{NULL} to disable.
#' @param ebWinsorize Logical. If TRUE, outlying first-pass raw values are
#' clamped to the outlier bounds before computing the EB prior. If FALSE, they
#' are set to missing for EB prior construction.
#' @param minsd Lower bound used for empirical raw SDs before model adjustment.
#' @param verbose Integer from 0 to 2. Passed to \code{\link{ctFit}}.
#' @param progress Logical. If TRUE, report the current EB fitting stage and
#' overwrite a single console line with the subject fitting percentage.
#' @param ... Additional arguments passed to each \code{\link{ctFit}} call.
#'
#' @return Object of class \code{ctEmpiricalBayesFit}, containing the subject
#' fit lists and metadata. \code{$initialfits} contains the first-pass model
#' prior fits, \code{$fits} contains the final empirical Bayes prior fits, and
#' \code{$passfits} contains every pass. Use \code{summary()} to compute final
#' raw-parameter means, SDs, covariance, and model summaries.
#' @export
#'
#' @examples
#' \donttest{
#' model <- ctModel(type='ct', manifestNames='Y1', latentNames='eta1',
#'   LAMBDA=matrix(1), silent=TRUE)
#' eb <- ctEmpiricalBayesFit(ctExample1, model,
#'   subjectFitArgs=list(optimcontrol=list(finishsamples=20)))
#' ebs <- summary(eb)
#' }
ctEmpiricalBayesFit <- function(datalong, ctstanmodel, subjects='all',
  priors=TRUE, optimize=TRUE, cores=1, subjectFitArgs=list(), Npasses=2,
  ebUse=c('rawest','rawposterior'), ebRobust=TRUE, ebOutlierMAD=6, ebOutlierQuantiles=c(.025,.975),
  ebWinsorize=TRUE, minsd=1e-6, verbose=0, progress=TRUE, ...){
  
  if(!'ctStanModel' %in% class(ctstanmodel)) stop('ctstanmodel must be a ctStanModel object')
  if(ctstanmodel$n.TIpred > 0 || length(ctstanmodel$TIpredNames) > 0) {
    stop('Time independent predictors are not supported for empirical Bayes subject-wise fitting')
  }
  
  datalong <- data.frame(datalong)
  if(!ctstanmodel$subjectIDname %in% colnames(datalong)) {
    stop('Subject id column ', ctstanmodel$subjectIDname, ' not found in data')
  }
  
  allsubjects <- unique(datalong[[ctstanmodel$subjectIDname]])
  if(length(subjects) == 1 && subjects %in% 'all') subjects <- allsubjects
  missingsubjects <- subjects[!subjects %in% allsubjects]
  if(length(missingsubjects) > 0) {
    stop('Requested subjects not found in data: ', paste(missingsubjects, collapse=', '))
  }
  if(length(subjects) < 2) stop('At least two subjects are required for empirical Bayes summaries')
  ebUse <- match.arg(ebUse)
  Npasses <- suppressWarnings(as.integer(Npasses[1]))
  if(!is.finite(Npasses) || is.na(Npasses) || Npasses < 2) {
    stop('Npasses must be an integer of at least 2')
  }
  
  subjectmodel <- ctstanmodel
  subjectmodel$pars$indvarying <- FALSE
  tieffects <- colnames(subjectmodel$pars)[grep('_effect', colnames(subjectmodel$pars), fixed=TRUE)]
  if(length(tieffects) > 0) subjectmodel$pars[, tieffects] <- FALSE
  
  dots <- list(...)
  fitargs <- list(
    ctstanmodel=subjectmodel,
    priors=priors,
    optimize=optimize,
    verbose=verbose)
  fitargs <- utils::modifyList(fitargs, subjectFitArgs)
  fitargs <- utils::modifyList(fitargs, dots)
  fitargs$cores <- 1
  fitargs <- ctEBfitArgsOptimDefaults(fitargs, stochastic=FALSE,
    firstpass=FALSE)
  initialfitargs <- ctEBfitArgsOptimDefaults(fitargs, stochastic=FALSE,
    firstpass=TRUE)
  
  initialfits <- ctEBfitSubjects(subjects=subjects, datalong=datalong,
    subjectIDname=ctstanmodel$subjectIDname, fitargs=initialfitargs,
    cores=cores, verbose=verbose, pass='model prior', progress=progress)
  
  parnames <- ctEBrawParnames(initialfits[[1]])
  rawlengths <- vapply(initialfits, function(fit) length(fit$stanfit$rawest), numeric(1))
  if(any(rawlengths != length(parnames))) {
    stop('Subject fits returned differing raw parameter counts')
  }
  
  passfits <- vector('list', Npasses)
  passmodels <- vector('list', Npasses)
  passraw <- vector('list', Npasses)
  passoriginalraw <- vector('list', Npasses)
  passrawForEB <- vector('list', Npasses - 1)
  passrawstats <- vector('list', Npasses - 1)
  passoutliers <- vector('list', Npasses - 1)
  passrawmaps <- vector('list', Npasses)
  names(passfits) <- names(passmodels) <- names(passraw) <-
    names(passoriginalraw) <- names(passrawmaps) <- paste0('pass', seq_len(Npasses))
  names(passrawForEB) <- names(passrawstats) <- names(passoutliers) <-
    paste0('pass', seq_len(Npasses - 1))
  
  passfits[[1]] <- initialfits
  passmodels[[1]] <- subjectmodel
  passrawmaps[[1]] <- ctEBidentityRawMap(parnames)
  currentfits <- initialfits
  currentrawmap <- passrawmaps[[1]]
  ebsubjectmodel <- adjustedmodel <- NULL
  
  for(passi in seq_len(Npasses - 1)){
    rawlocal <- ctEBrawMatrix(currentfits, parnames=parnames, use=ebUse)
    raworiginal <- ctEBmapRaw(rawlocal, currentrawmap)
    priorstats <- ctEBpriorRawStats(raworiginal, ebRobust=ebRobust,
      ebOutlierMAD=ebOutlierMAD, ebOutlierQuantiles=ebOutlierQuantiles,
      ebWinsorize=ebWinsorize, minsd=minsd, probs=c(.025,.5,.975))
    
    passraw[[passi]] <- rawlocal
    passoriginalraw[[passi]] <- raworiginal
    passrawForEB[[passi]] <- priorstats$rawForEB
    passrawstats[[passi]] <- priorstats$rawstats
    passoutliers[[passi]] <- priorstats$outliers
    
    ebsubjectmodel <- ctEBadjustModel(subjectmodel, priorstats$rawstats,
      sdscale='unit', minsd=minsd)
    adjustedmodel <- ebsubjectmodel
    
    ebfitargs <- fitargs
    ebfitargs$ctstanmodel <- ebsubjectmodel
    ebfitargs$inits <- NULL
    ebfitargs <- ctEBfitArgsOptimDefaults(ebfitargs, stochastic=FALSE,
      firstpass=FALSE)
    currentfits <- ctEBfitSubjects(subjects=subjects, datalong=datalong,
      subjectIDname=ctstanmodel$subjectIDname, fitargs=ebfitargs,
      cores=cores, verbose=verbose, pass=paste0('empirical Bayes prior ', passi),
      progress=progress)
    
    passfits[[passi + 1]] <- currentfits
    passmodels[[passi + 1]] <- ebsubjectmodel
    currentrawmap <- priorstats$rawstats[,c('param','mean','sd'),drop=FALSE]
    passrawmaps[[passi + 1]] <- currentrawmap
  }
  
  fits <- currentfits
  finalrawlocal <- ctEBrawMatrix(fits, parnames=parnames, use=ebUse)
  passraw[[Npasses]] <- finalrawlocal
  passoriginalraw[[Npasses]] <- ctEBmapRaw(finalrawlocal, currentrawmap)
  
  initialraw <- passoriginalraw[[1]]
  initialrawForEB <- passrawForEB[[1]]
  initialrawstats <- passrawstats[[1]]
  initialoutliers <- passoutliers[[1]]
  
  out <- list(
    call=match.call(),
    subjects=subjects,
    cores=cores,
    initialfits=initialfits,
    fits=fits,
    parnames=parnames,
    ebUse=ebUse,
    initialraw=initialraw,
    initialrawForEB=initialrawForEB,
    initialrawstats=initialrawstats,
    initialoutliers=initialoutliers,
    Npasses=Npasses,
    passfits=passfits,
    passmodels=passmodels,
    passraw=passraw,
    passoriginalraw=passoriginalraw,
    passrawForEB=passrawForEB,
    passrawstats=passrawstats,
    passoutliers=passoutliers,
    passrawmaps=passrawmaps,
    ebRobust=ebRobust,
    ebOutlierMAD=ebOutlierMAD,
    ebOutlierQuantiles=ebOutlierQuantiles,
    ebWinsorize=ebWinsorize,
    ctstanmodel=ctstanmodel,
    subjectmodel=subjectmodel,
    ebsubjectmodel=ebsubjectmodel,
    adjustedmodel=adjustedmodel,
    fitargs=fitargs,
    initialfitargs=initialfitargs,
    ebfitargs=ebfitargs)
  class(out) <- 'ctEmpiricalBayesFit'
  out
}

#' Summarise empirical Bayes subject-wise ctsem fits
#'
#' @param object Object returned by \code{\link{ctEmpiricalBayesFit}}.
#' @param use \code{'rawest'} to summarise final-pass subject point estimates,
#' or \code{'rawposterior'} to pool final-pass subject raw posterior samples.
#' @param probs Quantiles to report for raw parameters.
#' @param sdscale How to set \code{model$pars$sdscale} when reconstructing the
#' adjusted single-subject empirical Bayes model from the final empirical
#' raw distribution. \code{'unit'} keeps any later random-effect SDs on the
#' EB-standardised raw scale. \code{'rawsd'} uses the final empirical raw
#' SDs directly.
#' @param minsd Lower bound used for empirical raw SDs before model adjustment.
#' @param digits Number of digits for printed summary tables.
#' @param ... Unused.
#'
#' @return List containing final EB-prior raw estimates/samples, first-pass and
#' final original-raw parameter summaries, first-pass outlier diagnostics,
#' covariance/correlation matrices, the single-subject EB-adjusted model, and
#' the subject fit lists. \code{$raw} is on the final pass local raw scale;
#' \code{$originalraw} maps it back to the original model raw scale.
#' @method summary ctEmpiricalBayesFit
#' @export
summary.ctEmpiricalBayesFit <- function(object, use=c('rawest','rawposterior'),
  probs=c(.025,.5,.975), sdscale=c('unit','rawsd'), minsd=1e-6, digits=4, ...){
  
  if(!'ctEmpiricalBayesFit' %in% class(object)) {
    stop('Not a ctEmpiricalBayesFit object')
  }
  use <- match.arg(use)
  sdscale <- match.arg(sdscale)
  
  parnames <- object$parnames
  if(is.null(parnames)) parnames <- ctEBrawParnames(object$fits[[1]])
  initialraw <- object$initialraw
  if(is.null(initialraw)) initialraw <- ctEBrawMatrix(object$initialfits, parnames=parnames, use=object$ebUse)
  initialrawForEB <- object$initialrawForEB
  if(is.null(initialrawForEB)) initialrawForEB <- initialraw
  initialrawstats <- object$initialrawstats
  if(is.null(initialrawstats)) {
    initialrawstats <- ctEBrawStats(initialrawForEB, probs=probs)
    initialrawstats$sd[!is.finite(initialrawstats$sd) | initialrawstats$sd < minsd] <- minsd
  }
  passrawstats <- object$passrawstats
  if(is.null(passrawstats)) passrawstats <- list(initialrawstats)
  lastpriorstats <- passrawstats[[length(passrawstats)]]
  
  raw <- ctEBrawMatrix(object$fits, parnames=parnames, use=use)
  rawmap <- NULL
  if(!is.null(object$passrawmaps)) rawmap <- object$passrawmaps[[length(object$passrawmaps)]]
  originalraw <- ctEBmapRaw(raw, rawmap)
  rawstats <- ctEBrawStats(originalraw, probs=probs)
  rawstats$sd[!is.finite(rawstats$sd) | rawstats$sd < minsd] <- minsd
  
  adjustedmodel <- ctEBadjustModel(object$subjectmodel, lastpriorstats,
    sdscale=sdscale, minsd=minsd)
  
  initialrawcov <- stats::cov(initialraw, use='pairwise.complete.obs')
  initialrawcor <- suppressWarnings(stats::cor(initialraw, use='pairwise.complete.obs'))
  rawcov <- stats::cov(originalraw, use='pairwise.complete.obs')
  rawcor <- suppressWarnings(stats::cor(originalraw, use='pairwise.complete.obs'))
  
  out <- list(
    raw=raw,
    originalraw=originalraw,
    rawstats=rawstats,
    initialraw=initialraw,
    initialrawForEB=initialrawForEB,
    initialrawstats=initialrawstats,
    initialoutliers=object$initialoutliers,
    passraw=object$passraw,
    passoriginalraw=object$passoriginalraw,
    passrawForEB=object$passrawForEB,
    passrawstats=passrawstats,
    passoutliers=object$passoutliers,
    passrawmaps=object$passrawmaps,
    initialrawcov=initialrawcov,
    initialrawcor=initialrawcor,
    rawcov=rawcov,
    rawcor=rawcor,
    adjustedmodel=adjustedmodel,
    ebsubjectmodel=adjustedmodel,
    initialfits=object$initialfits,
    fits=object$fits,
    use=use,
    ebUse=object$ebUse,
    sdscale=sdscale,
    note='passrawstats contains the original-raw marginal summaries used to build each subsequent EB prior. raw is on the final pass local raw scale; originalraw maps it back to the original model raw scale. adjustedmodel is the random-effect-free single-subject model used for the final EB fitting pass. EB covariance/correlation are reported for inspection but only marginal raw means/SDs are represented in adjustedmodel.')
  class(out) <- 'summary.ctEmpiricalBayesFit'
  
  rounddf <- function(d){
    if(is.null(d)) return(NULL)
    data.frame(lapply(d, function(x){
      if(is.numeric(x)) round(x, digits) else x
    }), check.names=FALSE)
  }
  out$rawstats <- rounddf(out$rawstats)
  out$initialrawstats <- rounddf(out$initialrawstats)
  out$initialoutliers <- rounddf(out$initialoutliers)
  out$passrawstats <- lapply(out$passrawstats, rounddf)
  out$passoutliers <- lapply(out$passoutliers, rounddf)
  out$initialrawcov <- round(out$initialrawcov, digits)
  out$initialrawcor <- round(out$initialrawcor, digits)
  out$rawcov <- round(out$rawcov, digits)
  out$rawcor <- round(out$rawcor, digits)
  out
}
