# Function to compute the Hessian using bootstrap resampling
bootstrapHessian <- function(standata, sm, est, finishsamples, cores, scores=NULL) {
  # `scores` may be supplied by a backend that computes per-subject gradients
  # itself (see .ctBackendScoreMatrix); only the Stan path has to reconstruct
  # them a subject at a time.
  if(is.null(scores)) scores <- scorecalc(standata = standata,est = est,stanmodel = sm,
    subjectsonly = ctOptimNSubjects(standata) >= 2,
    returnsubjectlist = F,cores=cores)
  num_bootstrap_samples <- max(c(finishsamples,1000))
  alpha_max = 100 # Maximum bootstrap sample size factor
  alpha_min = 1 # Minimum bootstrap sample size factor
  n_threshold=1000 # Threshold n for alpha correction
  alpha <- alpha_max - (alpha_max - alpha_min) * (min(1000,nrow(scores))  / n_threshold)  # Bootstrap sample size factor
  num_bootstrap_samples  # Total number of bootstrap samples
  n <- nrow(scores)  # Number of observations
  p <- ncol(scores)  # Number of parameters
  
  # Create a bootstrap resampling matrix
  resample_matrix <- matrix(sample(1:n, size = round(alpha * n) * num_bootstrap_samples, replace = TRUE),
    nrow = num_bootstrap_samples, ncol = round(alpha * n))
  
  # Generate random weights for smoothing
  weights <- matrix(runif(length(resample_matrix), min = 0.1, max = 2),
    nrow = num_bootstrap_samples)
  
  # Aggregate gradients using matrix multiplication
  gradsamples <- matrix(0, nrow = num_bootstrap_samples, ncol = p)  # Initialize gradsamples
  
  # Compute the weighted sum of gradients for each bootstrap sample
  for (i in 1:num_bootstrap_samples) {
    gradsamples[i, ] <- colSums(scores[resample_matrix[i, ], , drop = FALSE] * weights[i, ])
  }
  
  # Trim outliers
  trim_percent <- 0#.05  # Proportion of outliers to trim
  if (trim_percent > 0) {
    lower <- apply(gradsamples, 2, quantile, probs = trim_percent, na.rm = TRUE)
    upper <- apply(gradsamples, 2, quantile, probs = 1 - trim_percent, na.rm = TRUE)
    gradsamples <- pmax(pmin(gradsamples, upper), lower)        }
  
  #  Compute the Hessian with alpha correction
  hess <- -corpcor::cov.shrink(gradsamples,verbose=FALSE) / alpha
  return(list(hess=hess,scores=scores))
}

ctOptimSafeCov <- function(cov, ridge=1e-8){
  cov <- as.matrix(cov)
  cov <- (cov + t(cov)) / 2
  if(any(!is.finite(cov))) stop('Non-finite covariance values')
  diagnostics <- list(nearPD=FALSE, ridgeApplied=FALSE,
    minEigenOriginal=NA_real_, minEigenFinal=NA_real_, ridge=ridge)
  eig <- try(eigen(cov, symmetric=TRUE), silent=TRUE)
  if('try-error' %in% class(eig)) {
    cov <- as.matrix(Matrix::nearPD(cov, conv.norm.type='F')$mat)
    eig <- eigen(cov, symmetric=TRUE)
    diagnostics$nearPD <- TRUE
  }
  mineig <- min(eig$values)
  diagnostics$minEigenOriginal <- mineig
  if(mineig <= ridge){
    eig$values <- pmax(eig$values, ridge)
    cov <- eig$vectors %*% diag(eig$values, length(eig$values)) %*%
      t(eig$vectors)
    cov <- (cov + t(cov)) / 2
    diagnostics$ridgeApplied <- TRUE
  }
  diagnostics$minEigenFinal <- min(eigen(cov, symmetric=TRUE,
    only.values=TRUE)$values)
  attr(cov, 'ctOptimSafeCov') <- diagnostics
  cov
}

ctOptimCovFromHessian <- function(hess, ridge=1e-8, warn=TRUE,
  context='Hessian'){
  hess <- (hess + t(hess)) / 2
  info <- -hess
  infoEig <- try(eigen(info, symmetric=TRUE, only.values=TRUE), silent=TRUE)
  minInfoEig <- if('try-error' %in% class(infoEig)) NA_real_ else
    min(infoEig$values)
  # The largest as well, because the smallest on its own says nothing. A
  # minimum of -7e-11 is rounding when the largest is 2.6e4 and a real rank
  # deficiency when the largest is 1e-9, and the repair warning below could not
  # tell those apart: it reported a magnitude with no scale to read it against.
  # A warning that fires identically on both is one people learn to ignore,
  # which is the worst outcome, because the text is the same when it mattered.
  maxInfoEig <- if('try-error' %in% class(infoEig)) NA_real_ else
    max(infoEig$values)
  infoEigenRatio <- if(is.finite(maxInfoEig) && maxInfoEig > 0 &&
      is.finite(minInfoEig)) abs(minInfoEig) / maxInfoEig else NA_real_
  covOk <- function(x){
    if('try-error' %in% class(x) || any(!is.finite(x))) return(FALSE)
    cholcheck <- try(suppressWarnings(chol((x + t(x)) / 2)), silent=TRUE)
    !'try-error' %in% class(cholcheck)
  }
  nearPDCov <- function(x){
    out <- try(Matrix::nearPD((x + t(x)) / 2, conv.norm.type='F',
        base.matrix=TRUE)$mat, silent=TRUE)
    if('try-error' %in% class(out)) return(out)
    (out + t(out)) / 2
  }
  repairSteps <- character()
  rawSolveSucceeded <- FALSE
  rawCholSucceeded <- FALSE
  usedNearPD <- FALSE
  usedInfoRidge <- FALSE
  usedGinv <- FALSE
  infoNearPD <- FALSE
  covNearPD <- FALSE
  covRidgeApplied <- FALSE
  covReady <- FALSE
  minInfoEigenFinal <- minInfoEig
  minCovEigenOriginal <- NA_real_
  minCovEigenFinal <- NA_real_
  
  rawcov <- try(suppressWarnings(solve(info)), silent=TRUE)
  rawSolveSucceeded <- !'try-error' %in% class(rawcov) &&
    all(is.finite(rawcov))
  if(rawSolveSucceeded) {
    rawcov <- (rawcov + t(rawcov)) / 2
    minCovEigenOriginal <- min(eigen(rawcov, symmetric=TRUE,
      only.values=TRUE)$values)
    rawCholSucceeded <- covOk(rawcov)
    if(rawCholSucceeded) {
      cov <- rawcov
      minCovEigenFinal <- minCovEigenOriginal
      diagnostics <- list(context=context, ridge=ridge,
        method='solve', minInfoEigenOriginal=minInfoEig,
        rawSolveSucceeded=rawSolveSucceeded,
        rawCholSucceeded=rawCholSucceeded,
        infoNearPD=FALSE, infoRidgeApplied=FALSE,
        minInfoEigenFinal=minInfoEigenFinal,
        usedNearPD=FALSE, usedGinv=FALSE,
        covNearPD=FALSE, covRidgeApplied=FALSE,
        minCovEigenOriginal=minCovEigenOriginal,
        minCovEigenFinal=minCovEigenFinal, repairSteps=repairSteps)
      attr(cov, 'ctOptimCovFromHessian') <- diagnostics
      return(cov)
    }
    repairSteps <- c(repairSteps,
      'solve(-hessian) succeeded but covariance was not positive definite')
    npdcov <- nearPDCov(rawcov)
    if(covOk(npdcov)) {
      cov <- npdcov
      covReady <- TRUE
      usedNearPD <- TRUE
      covNearPD <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps, 'nearPD applied to solved covariance')
    }
  } else {
    repairSteps <- c(repairSteps, 'solve(-hessian) failed')
  }
  
  if(!covReady) {
    safeInfo <- ctOptimSafeCov(info, ridge=ridge)
    infoDiagnostics <- attr(safeInfo, 'ctOptimSafeCov')
    infoNearPD <- isTRUE(infoDiagnostics$nearPD)
    usedInfoRidge <- isTRUE(infoDiagnostics$ridgeApplied)
    minInfoEigenFinal <- infoDiagnostics$minEigenFinal
    ridgecov <- try(suppressWarnings(solve(safeInfo)), silent=TRUE)
    if(covOk(ridgecov)) {
      cov <- (ridgecov + t(ridgecov)) / 2
      covReady <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps,
        'information matrix repaired before inversion')
    } else {
      repairSteps <- c(repairSteps,
        'solve() failed or covariance was not positive definite after information repair')
    }
  }
  
  if(!covReady) {
    usedGinv <- TRUE
    ginvcov <- try(MASS::ginv(info), silent=TRUE)
    if(!'try-error' %in% class(ginvcov)) ginvcov <- (ginvcov + t(ginvcov)) / 2
    if(covOk(ginvcov)) {
      cov <- ginvcov
      covReady <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps, 'MASS::ginv(-hessian) used')
    } else {
      npdcov <- if('try-error' %in% class(ginvcov)) ginvcov else
        nearPDCov(ginvcov)
      if(!covOk(npdcov)) {
        if('try-error' %in% class(ginvcov)) {
          stop('Could not construct covariance from Hessian using solve, ridge repair, or generalized inverse.',
            call.=FALSE)
        }
        npdcov <- ctOptimSafeCov(ginvcov, ridge=ridge)
        covRidgeApplied <- isTRUE(attr(npdcov, 'ctOptimSafeCov')$ridgeApplied)
      }
      cov <- npdcov
      covReady <- TRUE
      covNearPD <- TRUE
      minCovEigenFinal <- min(eigen(cov, symmetric=TRUE,
        only.values=TRUE)$values)
      repairSteps <- c(repairSteps,
        'MASS::ginv(-hessian) used with covariance positive-definite cleanup')
    }
  }
  
  if(is.na(minCovEigenOriginal) && covReady) {
    minCovEigenOriginal <- min(eigen((cov + t(cov)) / 2, symmetric=TRUE,
      only.values=TRUE)$values)
  }
  diagnostics <- list(context=context, ridge=ridge,
    minInfoEigenOriginal=minInfoEig,
    maxInfoEigenOriginal=maxInfoEig,
    infoEigenRatio=infoEigenRatio,
    # `sqrt(eps)` is where a symmetric eigendecomposition stops being able to
    # tell a small eigenvalue from zero, so below it the repair is arithmetic
    # rather than a statement about the model.
    #
    # The sign is half the test, and leaving it out was wrong. Rounding error
    # takes an eigenvalue that should be positive and makes it slightly
    # *negative*; it does not make it exactly zero. So a zero -- or a positive
    # value small enough to have needed the ridge -- is a genuinely singular
    # direction, which is the most serious case rather than the least, and the
    # ratio test alone classified it as negligible. `-diag(c(1, 0))` is the
    # minimal example, and it stopped warning.
    infoRepairNegligible=is.finite(infoEigenRatio) && minInfoEig < 0 &&
      infoEigenRatio < sqrt(.Machine$double.eps),
    rawSolveSucceeded=rawSolveSucceeded,
    rawCholSucceeded=rawCholSucceeded,
    infoNearPD=infoNearPD,
    infoRidgeApplied=usedInfoRidge,
    minInfoEigenFinal=minInfoEigenFinal,
    usedNearPD=usedNearPD,
    usedGinv=usedGinv,
    covNearPD=covNearPD,
    covRidgeApplied=covRidgeApplied,
    minCovEigenOriginal=minCovEigenOriginal,
    minCovEigenFinal=minCovEigenFinal,
    method=if(usedGinv) 'ginv' else if(usedInfoRidge) 'ridge_info'
      else if(usedNearPD) 'nearPD_cov' else 'solve',
    repairSteps=repairSteps)
  attr(cov, 'ctOptimCovFromHessian') <- diagnostics
  issues <- repairSteps
  if(isTRUE(diagnostics$infoNearPD)) issues <- c(issues,
    'nearPD was needed for the information matrix')
  if(isTRUE(diagnostics$infoRidgeApplied)) issues <- c(issues,
    paste0('information eigenvalues were floored at ridge=', ridge,
      ' (minimum original eigenvalue=', signif(minInfoEig, 4),
      if(is.finite(infoEigenRatio))
        paste0(', ', signif(infoEigenRatio, 3), ' of the largest)') else ')'))
  if(isTRUE(diagnostics$usedGinv)) issues <- c(issues,
    'MASS::ginv() was used')
  if(isTRUE(diagnostics$covNearPD) || isTRUE(diagnostics$covRidgeApplied)) {
    issues <- c(issues,
      'the resulting covariance required positive-definite cleanup')
  }
  if(warn && length(issues) > 0) {
    # Graded. A repair that floored an eigenvalue indistinguishable from zero is
    # arithmetic, and is reported as such; one that floored a substantively
    # negative or tiny eigenvalue is a statement about what the data can
    # determine, and keeps the warning. Making that distinction is the point --
    # the same text for both is what taught people to ignore it, and the cases
    # it currently conflates are genuinely different. A near-integrated trend
    # process *should* warn here.
    if(isTRUE(diagnostics$infoRepairNegligible) &&
        !isTRUE(diagnostics$usedGinv) && !isTRUE(diagnostics$infoNearPD)) {
      # Deliberately not "arithmetic, not a statement about the model", which
      # is what this said and could not support. A smallest eigenvalue of
      # -1e-12 relative to the largest is what rounding does to a
      # positive-definite matrix, and it is *also* what rounding does to a
      # genuinely singular one -- the two are indistinguishable from this
      # number alone. Usually the first, so a message rather than a warning;
      # never certainly the first, so it says which is which is checkable
      # elsewhere rather than pronouncing.
      message(context, ' covariance: the information matrix needed a numerical ',
        'nudge before inversion (smallest eigenvalue ',
        signif(minInfoEig, 3), ', ', signif(infoEigenRatio, 3),
        ' of the largest, which is indistinguishable from zero at machine ',
        'precision). Usually rounding rather than a flat direction, but the ',
        'two look the same at this magnitude; fit$identifiability names the ',
        'parameters if it is the latter.')
    } else {
      warning(context, ' covariance from Hessian required numerical repair: ',
        paste(issues, collapse='; '),
        if(is.finite(infoEigenRatio) &&
            infoEigenRatio >= sqrt(.Machine$double.eps))
          paste0('. The floored eigenvalue is ', signif(infoEigenRatio, 3),
            ' of the largest, too large to be rounding: some direction of this ',
            'model is close to unidentified and the standard errors along it ',
            'are not trustworthy') else '',
        call.=FALSE)
    }
  }
  cov
}

ctOptimNormalDraws <- function(mean, cov, n, df=Inf){
  cov <- ctOptimSafeCov(cov)
  z <- matrix(stats::rnorm(n * length(mean)), nrow=n)
  draws <- z %*% chol(cov)
  if(is.finite(df)){
    draws <- draws / sqrt(stats::rchisq(n, df=df) / df)
  }
  sweep(draws, 2, mean, '+')
}

ctOptimScoreMatrix <- function(standata, sm, est, cores=1, scores=NULL){
  if(!is.null(scores)) return(scores)
  scorecalc(standata=standata, est=est, stanmodel=sm,
    subjectsonly=ctOptimNSubjects(standata) >= 2,
    returnsubjectlist=FALSE,
    cores=cores)
}

ctOptimNSubjects <- function(standata){
  nsubjects <- suppressWarnings(as.integer(standata$nsubjects[1]))
  if(length(nsubjects) < 1 || is.na(nsubjects) || !is.finite(nsubjects)) {
    if(!is.null(standata$subject)) nsubjects <- length(unique(standata$subject))
  }
  if(length(nsubjects) < 1 || is.na(nsubjects) || !is.finite(nsubjects)) {
    nsubjects <- NA_integer_
  }
  nsubjects
}

ctOptimNDataPoints <- function(standata){
  ndatapoints <- suppressWarnings(as.integer(standata$ndatapoints[1]))
  if(length(ndatapoints) < 1 || is.na(ndatapoints) || !is.finite(ndatapoints)) {
    if(!is.null(standata$subject)) ndatapoints <- length(standata$subject)
  }
  if(length(ndatapoints) < 1 || is.na(ndatapoints) || !is.finite(ndatapoints)) {
    ndatapoints <- NA_integer_
  }
  ndatapoints
}

ctOptimCheckUncertaintyData <- function(standata, uncertainty, finishsamples,
  npars=NULL){
  nsubjects <- ctOptimNSubjects(standata)
  ndatapoints <- ctOptimNDataPoints(standata)
  if(is.null(npars)) npars <- NA_integer_
  npars <- suppressWarnings(as.integer(npars[1]))
  if(length(npars) < 1 || is.na(npars) || !is.finite(npars)) npars <- NA_integer_
  if(uncertainty %in% c('bootstrap','fullbootstrap') &&
      finishsamples < 2) {
    stop(uncertainty, ' uncertainty requires at least two samples / refits ',
      'to estimate a covariance; increase finishsamples.', call.=FALSE)
  }
  if(uncertainty == 'fullbootstrap'){
    if(is.na(nsubjects) || nsubjects < 2) {
      stop('fullbootstrap uncertainty requires at least two subjects.',
        call.=FALSE)
    }
    if(nsubjects < 10) {
      warning('fullbootstrap uncertainty requested with fewer than ten ',
        'independent subjects; the bootstrap distribution may be unstable.',
        call.=FALSE)
    }
    if(!is.na(npars) && finishsamples <= npars) {
      warning('fullbootstrap requested with finishsamples <= number of raw ',
        'parameters; the empirical covariance is rank limited and will be ',
        'regularised.', call.=FALSE)
    }
  }
  if(uncertainty %in% c('bootstrap','sandwich','opg')){
    subjectScores <- !is.na(nsubjects) && nsubjects >= 2
    nscore <- if(subjectScores) nsubjects else ndatapoints
    if(is.na(nscore) || nscore < 2) {
      stop(uncertainty, ' uncertainty requires at least two ',
        if(subjectScores) 'subject-level' else 'case-level',
        ' score contribution rows.', call.=FALSE)
    }
    if(subjectScores && nscore < 10) {
      warning(uncertainty, ' uncertainty is based on fewer than ten ',
        'independent subject-level score contributions; the covariance ',
        'estimate may be unstable.', call.=FALSE)
    }
    if(!is.na(npars) && nscore <= npars) {
      warning(uncertainty, ' uncertainty has no more score contribution rows ',
        'than raw parameters; the covariance estimate is rank limited and ',
        'will be regularised.', call.=FALSE)
    }
    if(!subjectScores) {
      warning(uncertainty, ' uncertainty for a single-subject model uses ',
        'case-level score contributions. This can be unreliable when ',
        'observations are serially dependent.',
        call.=FALSE)
    }
  }
  invisible(list(nsubjects=nsubjects, ndatapoints=ndatapoints))
}

ctOptimBootstrapDraws <- function(est, cov, scores, n=1000){
  scores <- as.matrix(scores)
  scores <- scale(scores, center=TRUE, scale=FALSE)
  nscore <- nrow(scores)
  if(nscore < 2) stop('At least two score rows are required for bootstrap uncertainty')
  draws <- matrix(NA_real_, nrow=n, ncol=length(est))
  for(i in seq_len(n)){
    idx <- sample.int(nscore, nscore, replace=TRUE)
    score_sum <- colSums(scores[idx,,drop=FALSE])
    draws[i,] <- est + as.numeric(cov %*% score_sum)
  }
  draws
}

ctOptimBootstrapStandata <- function(standata, subjects){
  subjects <- as.integer(subjects)
  if(length(subjects) < 2) stop('At least two subjects are required for full bootstrap uncertainty')
  long <- standatatolong(standata)
  longlist <- vector('list', length(subjects))
  for(i in seq_along(subjects)){
    longi <- long[long$subject %in% subjects[i], , drop=FALSE]
    if(nrow(longi) < 1) stop('Subject ', subjects[i], ' not found in standata')
    longi$subject <- i
    longlist[[i]] <- longi
  }
  longboot <- do.call(rbind, longlist)
  row.names(longboot) <- NULL
  standataboot <- standatalongremerge(long=longboot, standata=standata)
  standataboot$ndatapoints <- as.integer(nrow(longboot))
  standataboot$nsubjects <- as.integer(length(subjects))
  standataboot$subject <- array(as.integer(longboot$subject))
  if(standata$ntipred > 0) {
    standataboot$tipredsdata <- standata$tipredsdata[subjects, , drop=FALSE]
  }
  standataboot$idmap <- data.frame(
    original=paste0('boot', seq_along(subjects), '_subject', subjects),
    new=seq_along(subjects))
  standataboot
}

ctOptimFullBootstrapOne <- function(i, est, standata, sm, fitCores, tol,
  verbose=0){
  subjects <- unique(standata$subject)
  sampledSubjects <- sample(subjects, length(subjects), replace=TRUE)
  standataboot <- ctOptimBootstrapStandata(standata=standata,
    subjects=sampledSubjects)
  standataboot$savesubjectmatrices <- 0L
  standataboot$nsubsets <- 1L
  lpgsetup <- ctOptimDataLpgFunc(sm=sm, standata=standataboot,
    cores=fitCores)
  on.exit({
    if(!is.null(lpgsetup$cl)) try(parallel::stopCluster(lpgsetup$cl),
      silent=TRUE)
    if(!is.null(lpgsetup$smfile) && nzchar(lpgsetup$smfile)) {
      try(file.remove(lpgsetup$smfile), silent=TRUE)
    }
  }, add=TRUE)
  if(verbose > 0) message('Full bootstrap sample ', i)
  opt <- try(ctOptim(init=est, lpgFunc=lpgsetup$lpg, tol=tol,
    nsubsets=1L, stochastic=FALSE, stochasticTolAdjust=1,
    bfgsType='mize'), silent=TRUE)
  if('try-error' %in% class(opt) || is.null(opt$par) ||
      length(opt$par) != length(est) || any(!is.finite(opt$par))) {
    msg <- if('try-error' %in% class(opt) && !is.null(attr(opt,
          'condition'))) {
      conditionMessage(attr(opt, 'condition'))
    } else 'non-finite optimized parameters'
    return(list(ok=FALSE, par=rep(NA_real_, length(est)),
      sampledSubjects=sampledSubjects, message=msg))
  }
  list(ok=TRUE, par=opt$par, sampledSubjects=sampledSubjects,
    value=opt$value, message=opt$message)
}

ctOptimFullBootstrapDraws <- function(est, standata, sm, n=1000, cores=1,
  control=list(), verbose=0){
  if(standata$nsubjects < 2) {
    stop('Full bootstrap uncertainty requires at least two subjects')
  }
  if(is.null(control$bootstrapFitCores)) control$bootstrapFitCores <- 1L
  if(is.null(control$bootstrapTol)) control$bootstrapTol <- 1e-5
  fitCores <- suppressWarnings(as.integer(control$bootstrapFitCores[1]))
  if(!is.finite(fitCores) || is.na(fitCores) || fitCores < 1) fitCores <- 1L
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  outerCores <- min(n, max(1L, floor(cores / fitCores)))
  fitCores <- min(fitCores, standata$nsubjects)
  if(verbose > 0 || outerCores > 1) {
    message('Fitting ', n, ' full bootstrap samples using ', outerCores,
      ' bootstrap worker(s) and ', fitCores, ' core(s) per refit')
  }
  
  if(outerCores > 1){
    cl <- makeClusterID(outerCores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
    bootHelpers <- c('ctOptimFullBootstrapOne',
      'ctOptimBootstrapStandata', 'ctOptimDataLpgFunc', 'ctOptim',
      'standatatolong', 'standatalongremerge', 'standatalongobjects',
      'stan_reinitsf', 'getcxxfun', 'suppressOutput', 'makeClusterID',
      'parallelStanSetup', 'clusterIDexport', 'clusterIDeval',
      'singlecoreStanSetup', 'parlptext')
    parallel::clusterExport(cl, c('bootHelpers', bootHelpers),
      envir=environment(ctOptimFullBootstrapDraws))
    parallel::clusterEvalQ(cl, {
      bootEnv <- new.env(parent=.GlobalEnv)
      for(fn in bootHelpers){
        obj <- get(fn, envir=.GlobalEnv)
        if(is.function(obj)){
          environment(obj) <- bootEnv
        }
        assign(fn, obj, envir=bootEnv)
      }
      options(ctsem.bootstrap.env=bootEnv)
      rm(list=c(bootHelpers, 'bootHelpers'), envir=.GlobalEnv)
      rm(bootHelpers, bootEnv)
      NULL
    })
    out <- parallel::parLapplyLB(cl, seq_len(n), function(i){
      bootEnv <- getOption('ctsem.bootstrap.env')
      if(is.null(bootEnv)) stop('Missing ctsem bootstrap worker environment')
      get('ctOptimFullBootstrapOne', envir=bootEnv)(i=i, est=est, standata=standata,
        sm=sm, fitCores=fitCores, tol=control$bootstrapTol,
        verbose=verbose)
    })
  } else {
    out <- lapply(seq_len(n), function(i){
      if(verbose == 0) message('\rFull bootstrap sample ', i, '/', n,
        appendLF=FALSE)
      ctOptimFullBootstrapOne(i=i, est=est, standata=standata, sm=sm,
        fitCores=fitCores, tol=control$bootstrapTol, verbose=verbose)
    })
    if(verbose == 0) message('')
  }
  ok <- vapply(out, `[[`, logical(1), 'ok')
  if(!any(ok)) {
    msgs <- unique(vapply(out, function(x) x$message, character(1)))
    stop('All full bootstrap refits failed. First errors: ',
      paste(utils::head(msgs, 3), collapse='; '))
  }
  if(any(!ok)) {
    warning(sum(!ok), ' full bootstrap refits failed and were omitted.',
      call.=FALSE)
  }
  draws <- do.call(rbind, lapply(out[ok], `[[`, 'par'))
  list(draws=draws,
    sampledSubjects=lapply(out[ok], `[[`, 'sampledSubjects'),
    failures=out[!ok], outerCores=outerCores, fitCores=fitCores,
    tol=control$bootstrapTol)
}

ctOptimSurrogateDirections <- function(p, n){
  dirs <- rbind(diag(p), -diag(p))
  while(nrow(dirs) < n){
    addn <- ceiling((n - nrow(dirs)) / 2)
    z <- matrix(stats::rnorm(addn * p), nrow=addn)
    z <- z / sqrt(rowSums(z^2))
    dirs <- rbind(dirs, z, -z)
  }
  dirs[seq_len(n),,drop=FALSE]
}

ctOptimSurrogateEvalPoint <- function(est, lpgFunc, cholcov, z, baseValue){
  rawstep <- as.numeric(z %*% cholcov)
  lp <- try(suppressMessages(suppressWarnings(lpgFunc(est + rawstep))),
    silent=TRUE)
  value <- NA_real_
  grad <- rep(NA_real_, length(est))
  if(!'try-error' %in% class(lp)) {
    value <- lp[1]
    lpgrad <- attributes(lp)$gradient
    if(length(lpgrad) == length(est)) grad <- lpgrad
  }
  drop <- baseValue - value
  finite <- is.finite(value) && is.finite(drop) && all(is.finite(grad))
  list(white=as.numeric(z), raw=rawstep, value=value, gradient=grad,
    drop=drop, finite=finite)
}

ctOptimSurrogateBestEval <- function(evals, targetDrop, preferRange=NULL){
  drops <- vapply(evals, `[[`, numeric(1), 'drop')
  finite <- vapply(evals, `[[`, logical(1), 'finite') &
    is.finite(drops) & drops > 0
  if(!is.null(preferRange)) {
    inrange <- finite & drops >= preferRange[1] & drops <= preferRange[2]
    if(any(inrange)) finite <- inrange
  }
  if(!any(finite)) {
    finite <- vapply(evals, `[[`, logical(1), 'finite')
  }
  if(!any(finite)) return(evals[[length(evals)]])
  ii <- which(finite)
  ii <- ii[which.min(abs(log(pmax(drops[ii], .Machine$double.eps) /
      targetDrop)))]
  evals[[ii]]
}

ctOptimSurrogateRadiusForDrop <- function(radius, drop, targetDrop,
  minFactor=.2, maxFactor=8, safety=1.1){
  if(!is.finite(drop) || drop <= 0) return(radius * maxFactor)
  factor <- sqrt(targetDrop / drop) * safety
  radius * min(maxFactor, max(minFactor, factor))
}

ctOptimSurrogateTargetPoint <- function(est, lpgFunc, cholcov, direction,
  targetDrop, dropRange, initialRadius, maxRadius=64, maxEval=6,
  targetFactor=1.35, baseValue){
  direction <- as.numeric(direction)
  direction <- direction / sqrt(sum(direction^2))
  evals <- list()
  evalAt <- function(radius){
    ev <- ctOptimSurrogateEvalPoint(est=est, lpgFunc=lpgFunc,
      cholcov=cholcov, z=radius * direction, baseValue=baseValue)
    ev$radius <- radius
    ev
  }
  addEval <- function(radius){
    evals[[length(evals) + 1L]] <<- evalAt(radius)
    evals[[length(evals)]]
  }
  goodTarget <- function(ev){
    ev$finite && is.finite(ev$drop) && ev$drop > 0 &&
      abs(log(ev$drop / targetDrop)) <= log(targetFactor)
  }
  initialRadius <- max(.Machine$double.eps, initialRadius)
  ev <- addEval(initialRadius)
  if(goodTarget(ev)) {
    ev$targeted <- TRUE
    ev$neval <- length(evals)
    return(ev)
  }
  
  low <- 0
  high <- initialRadius
  if(ev$finite && is.finite(ev$drop) && ev$drop > 0 &&
      ev$drop < targetDrop) {
    low <- initialRadius
    while(length(evals) < maxEval && high < maxRadius) {
      high <- min(maxRadius, ctOptimSurrogateRadiusForDrop(high,
        ev$drop, targetDrop, minFactor=1.5))
      ev <- addEval(high)
      if(goodTarget(ev)) {
        ev$targeted <- TRUE
        ev$neval <- length(evals)
        return(ev)
      }
      if(!ev$finite || !is.finite(ev$drop) || ev$drop >= targetDrop) break
      low <- high
    }
  }
  
  while(length(evals) < maxEval && high > 0 && high > low) {
    if(is.finite(ev$drop) && ev$drop > targetDrop && low == 0) {
      mid <- ctOptimSurrogateRadiusForDrop(high, ev$drop, targetDrop,
        minFactor=.1, maxFactor=.8, safety=.9)
      mid <- min(high * .95, max(.Machine$double.eps, mid))
    } else {
      mid <- (low + high) / 2
    }
    ev <- addEval(mid)
    if(goodTarget(ev)) {
      ev$targeted <- TRUE
      ev$neval <- length(evals)
      return(ev)
    }
    if(!ev$finite || !is.finite(ev$drop) || ev$drop >= targetDrop) {
      high <- mid
    } else {
      low <- mid
    }
  }
  
  ev <- ctOptimSurrogateBestEval(evals=evals, targetDrop=targetDrop,
    preferRange=dropRange)
  ev$targeted <- ev$finite && is.finite(ev$drop) && ev$drop >= dropRange[1] &&
    ev$drop <= dropRange[2]
  ev$neval <- length(evals)
  ev
}

ctOptimSurrogateBacktransformHessian <- function(hessWhite, cholcov){
  invchol <- backsolve(cholcov, diag(ncol(cholcov)))
  hess <- invchol %*% hessWhite %*% t(invchol)
  (hess + t(hess)) / 2
}

ctOptimSurrogateProfileDirections <- function(est, lpgFunc, cholcov,
  directions, targetDrop=2, maxStep=64, tol=.02, maxIter=25,
  initialStep=1, maxExpand=4, baseValue=NULL, verbose=0){
  if(is.null(baseValue)) {
    baseValue <- suppressMessages(suppressWarnings(lpgFunc(est)))[1]
  }
  directions <- as.matrix(directions)
  if(nrow(directions) < 1) {
    return(data.frame())
  }
  evalDrop <- function(z){
    rawstep <- as.numeric(z %*% cholcov)
    lp <- try(suppressMessages(suppressWarnings(lpgFunc(est + rawstep))),
      silent=TRUE)
    if('try-error' %in% class(lp)) return(NA_real_)
    baseValue - lp[1]
  }
  out <- vector('list', nrow(directions) * 2L)
  oi <- 0L
  for(di in seq_len(nrow(directions))){
    diri <- directions[di,]
    diri <- diri / sqrt(sum(diri^2))
    maxStepi <- maxStep[min(length(maxStep), di)]
    initialStepi <- initialStep[min(length(initialStep), di)]
    if(!is.finite(maxStepi) || maxStepi <= 0) maxStepi <- 64
    if(!is.finite(initialStepi) || initialStepi <= 0) initialStepi <- 1
    initialStepi <- min(initialStepi, maxStepi)
    for(sgn in c(-1, 1)){
      oi <- oi + 1L
      low <- 0
      high <- initialStepi
      dhigh <- evalDrop(sgn * high * diri)
      while(is.finite(dhigh) && dhigh < targetDrop && high < maxStepi){
        low <- high
        high <- min(maxStepi, ctOptimSurrogateRadiusForDrop(high,
          dhigh, targetDrop, minFactor=1.5))
        dhigh <- evalDrop(sgn * high * diri)
      }
      expand <- 0L
      while(is.finite(dhigh) && dhigh < targetDrop && expand < maxExpand) {
        expand <- expand + 1L
        low <- high
        high2 <- ctOptimSurrogateRadiusForDrop(high, dhigh, targetDrop,
          minFactor=1.5, maxFactor=16)
        if(!is.finite(high2) || high2 <= high) high2 <- high * 2
        high <- high2
        dhigh <- evalDrop(sgn * high * diri)
      }
      reached <- is.finite(dhigh) && dhigh >= targetDrop
      if(reached) {
        for(iter in seq_len(maxIter)){
          if(low == 0 && is.finite(dhigh) && dhigh > targetDrop) {
            mid <- ctOptimSurrogateRadiusForDrop(high, dhigh, targetDrop,
              minFactor=.1, maxFactor=.8, safety=.9)
            mid <- min(high * .95, max(.Machine$double.eps, mid))
          } else {
            mid <- (low + high) / 2
          }
          dmid <- evalDrop(sgn * mid * diri)
          if(!is.finite(dmid) || dmid >= targetDrop) high <- mid else low <- mid
          if(is.finite(dmid) && abs(dmid - targetDrop) <= tol * targetDrop) break
        }
        step <- high
        drop <- evalDrop(sgn * step * diri)
        curvature <- 2 * targetDrop / step^2
      } else {
        step <- high
        drop <- dhigh
        curvature <- if(is.finite(drop) && drop > 0) 2 * drop / step^2
          else NA_real_
      }
      status <- if(reached) 'reached' else if(!is.finite(drop)) {
        'nonfinite'
      } else if(drop < targetDrop) {
        'belowTargetAtMaxStep'
      } else {
        'notReached'
      }
      out[[oi]] <- data.frame(direction=di, sign=sgn, step=step,
        drop=drop, curvature=curvature, reached=reached, status=status,
        expansions=expand)
    }
  }
  profiles <- do.call(rbind, out)
  if(verbose > 0) {
    message('Surrogate profiled ', nrow(directions), ' direction(s); ',
      sum(profiles$reached), '/', nrow(profiles),
      ' one-sided targets reached')
    if(any(!profiles$reached)) {
      message('Surrogate profile misses: ',
        paste(names(table(profiles$status[!profiles$reached])),
          as.integer(table(profiles$status[!profiles$reached])),
          sep='=', collapse=', '))
    }
  }
  profiles
}

ctOptimSurrogateProfileCurvature <- function(hessWhite, est, lpgFunc,
  cholcov, targetDrop=2, maxStep=64, baseValue=NULL, verbose=0){
  infoWhite <- -((hessWhite + t(hessWhite)) / 2)
  eig <- eigen(infoWhite, symmetric=TRUE)
  directions <- t(eig$vectors)
  expectedStep <- rep(maxStep, length(eig$values))
  positive <- is.finite(eig$values) & eig$values > 0
  expectedStep[positive] <- sqrt(2 * targetDrop / eig$values[positive])
  expectedStep[!is.finite(expectedStep) | expectedStep <= 0] <- maxStep
  profileMaxStep <- pmax(maxStep, expectedStep * 2)
  profileInitialStep <- pmax(.Machine$double.eps, expectedStep)
  profiles <- ctOptimSurrogateProfileDirections(est=est, lpgFunc=lpgFunc,
    cholcov=cholcov, directions=directions, targetDrop=targetDrop,
    maxStep=profileMaxStep, initialStep=profileInitialStep,
    baseValue=baseValue, verbose=verbose)
  adjusted <- 0L
  newvals <- eig$values
  for(i in seq_along(newvals)){
    curv <- profiles$curvature[profiles$direction == i &
        is.finite(profiles$curvature)]
    if(length(curv) > 0) {
      profileCurv <- max(curv)
      if(is.finite(profileCurv) && profileCurv > newvals[i]) {
        newvals[i] <- profileCurv
        adjusted <- adjusted + 1L
      }
    }
  }
  infoWhite <- eig$vectors %*% diag(newvals, length(newvals)) %*%
    t(eig$vectors)
  infoWhite <- (infoWhite + t(infoWhite)) / 2
  list(hessWhite=-infoWhite, profiles=profiles, nProfiled=length(newvals),
    nAdjusted=adjusted)
}

ctOptimSurrogateHessian <- function(est, lpgFunc, cov, npoints=NULL,
  scale=.5, ridge=1e-6, profile=TRUE, profileTargetDrop=NULL,
  profileMaxStep=64, verbose=0){
  p <- length(est)
  if(is.null(npoints)) npoints <- max(4 * p, 50)
  cov <- ctOptimSafeCov(cov, ridge=ridge)
  cholcov <- chol(cov)
  targetDrop <- 2
  dropRange <- c(.25, 6)
  defaultScale <- .5
  initialRadius <- (scale / defaultScale) * sqrt(2 * targetDrop)
  base <- suppressMessages(suppressWarnings(lpgFunc(est)))
  baseValue <- base[1]
  basegrad <- attributes(base)$gradient
  if(is.null(basegrad) || any(!is.finite(basegrad))) basegrad <- rep(0, p)
  
  directions <- ctOptimSurrogateDirections(p=p, n=npoints)
  evals <- vector('list', npoints)
  for(i in seq_len(npoints)){
    if(verbose > 0) message('\rFitting local quadratic surrogate, point ',
      i, '/', npoints, appendLF=FALSE)
    evals[[i]] <- ctOptimSurrogateTargetPoint(est=est, lpgFunc=lpgFunc,
      cholcov=cholcov, direction=directions[i,], targetDrop=targetDrop,
      dropRange=dropRange, initialRadius=initialRadius,
      baseValue=baseValue)
  }
  if(verbose > 0) message('')
  values <- vapply(evals, `[[`, numeric(1), 'value')
  drops <- vapply(evals, `[[`, numeric(1), 'drop')
  targeted <- vapply(evals, `[[`, logical(1), 'targeted')
  neval <- vapply(evals, `[[`, numeric(1), 'neval')
  gradients <- do.call(rbind, lapply(evals, `[[`, 'gradient'))
  design <- do.call(rbind, lapply(evals, `[[`, 'raw'))
  whiteDesign <- do.call(rbind, lapply(evals, `[[`, 'white'))
  finite <- is.finite(values) & is.finite(drops) &
    apply(gradients, 1, function(x) all(is.finite(x)))
  inrange <- finite & drops >= dropRange[1] & drops <= dropRange[2]
  positive <- finite & drops > 0
  keep <- inrange
  if(sum(keep) < npoints && any(positive & !keep)) {
    add <- which(positive & !keep)
    add <- add[order(abs(log(pmax(drops[add], .Machine$double.eps) /
        targetDrop)))]
    keep[add[seq_len(min(length(add), npoints - sum(keep)))]] <- TRUE
  }
  if(sum(keep) <= p) stop('Too few finite surrogate evaluations')
  design <- design[keep,,drop=FALSE]
  whiteDesign <- whiteDesign[keep,,drop=FALSE]
  gradients <- gradients[keep,,drop=FALSE]
  values <- values[keep]
  drops <- drops[keep]
  diagnostics <- list(nRequested=npoints, nUsed=nrow(design),
    nTargeted=sum(targeted), nInRange=sum(inrange),
    nFinitePositive=sum(positive), nFinite=sum(finite),
    meanEvaluations=mean(neval), maxEvaluations=max(neval),
    targetDrop=targetDrop, dropRange=dropRange)
  rawgrad <- sweep(gradients, 2, basegrad, '-')
  y <- rawgrad %*% t(cholcov)
  pointWeights <- 1 / pmax(abs(log(pmax(drops, .Machine$double.eps) /
        targetDrop)), .25)
  pointWeights <- pointWeights / mean(pointWeights, na.rm=TRUE)
  sqrtw <- sqrt(pointWeights)
  xw <- whiteDesign * sqrtw
  yw <- y * sqrtw
  xtx <- crossprod(xw) + diag(ridge, p)
  coef <- solve(xtx, crossprod(xw, yw))
  hessWhite <- (t(coef) + coef) / 2
  if(isTRUE(profile)) {
    profiled <- ctOptimSurrogateProfileCurvature(hessWhite=hessWhite,
      est=est, lpgFunc=lpgFunc, cholcov=cholcov,
      targetDrop=if(is.null(profileTargetDrop)) targetDrop else profileTargetDrop,
      maxStep=profileMaxStep, baseValue=baseValue, verbose=verbose)
    hessWhite <- profiled$hessWhite
  } else {
    profiled <- list(profiles=data.frame(), nProfiled=0L, nAdjusted=0L,
      note='Directional profiling disabled')
  }
  hess <- ctOptimSurrogateBacktransformHessian(hessWhite, cholcov)
  list(hessian=hess, values=values, gradients=gradients, design=design,
    whiteDesign=whiteDesign, drops=drops, targetDrop=targetDrop, dropRange=dropRange,
    scale=initialRadius, parScale=rep(1, p), nfinite=sum(keep),
    rounds=1, diagnostics=diagnostics, profile=profiled)
}

ctOptimComputeUncertainty <- function(est, standata, sm, lpgFunc,
  uncertainty=c('hessian','surrogate','is','bootstrap','fullbootstrap',
    'sandwich','opg'),
  finishsamples=1000, cores=1, matsetup=NA, control=list(), verbose=0,
  scores=NULL, hessian=NULL){
  
  uncertainty <- match.arg(uncertainty)
  ctOptimCheckUncertaintyData(standata=standata, uncertainty=uncertainty,
    finishsamples=finishsamples, npars=length(est))
  if(is.null(control$ridge)) control$ridge <- 1e-8
  if(is.null(control$hessianStep)) control$hessianStep <- 1e-3
  if(is.null(control$surrogateScale)) control$surrogateScale <- .5
  if(is.null(control$surrogateNpoints)) control$surrogateNpoints <- NULL
  if(is.null(control$surrogateProfile)) control$surrogateProfile <- TRUE
  if(is.null(control$surrogateProfileTargetDrop)) {
    control$surrogateProfileTargetDrop <- NULL
  }
  if(is.null(control$surrogateProfileMaxStep)) {
    control$surrogateProfileMaxStep <- 64
  }
  
  base <- suppressMessages(suppressWarnings(lpgFunc(est)))
  base_gradient <- attributes(base)$gradient
  if(is.null(base_gradient) || any(!is.finite(base_gradient))) {
    base_gradient <- rep(0, length(est))
  }
  
  hessian_result <- NULL
  scoremat <- NULL
  draws <- NULL
  method_details <- list()
  covavailable <- FALSE
  
  if(uncertainty == 'surrogate' && !is.null(control$initialCov)){
    cov <- ctOptimSafeCov(control$initialCov, ridge=control$ridge)
    covavailable <- TRUE
  }
  
  if(uncertainty %in% c('hessian','sandwich','bootstrap','is') ||
      (uncertainty == 'surrogate' && !covavailable)){
    if(!is.null(hessian)){
      # An exact Hessian supplied by the caller -- see .ctBackendHessian(),
      # which gets one from the julia engine by differentiating its own
      # reverse-mode gradient. There is nothing to difference, so there is no
      # step to choose and no pair of one-sided estimates to reconcile.
      hess <- hessian
    } else {
      message('Estimating Hessian')
      hess1 <- numericHessianFunc(pars=est, step=control$hessianStep,
        verbose=verbose, directions=1, lpgFunc=lpgFunc,
        base_value=base[1], base_gradient=base_gradient)
      hess2 <- numericHessianFunc(pars=est, step=control$hessianStep,
        verbose=verbose, directions=-1, lpgFunc=lpgFunc,
        base_value=base[1], base_gradient=base_gradient)
      message('')
      hessian_result <- processHessianMatrices(hess1, hess2, verbose, matsetup)
      hess <- hessian_result$hess
    }
    cov <- ctOptimCovFromHessian(hess, ridge=control$ridge)
    covavailable <- TRUE
  }
  
  if(uncertainty == 'opg'){
    message('Estimating score / OPG covariance')
    score_hessian <- bootstrapHessian(standata=standata, sm=sm, est=est,
      finishsamples=finishsamples, cores=cores, scores=scores)
    hess <- score_hessian$hess
    scoremat <- score_hessian$scores
    cov <- ctOptimCovFromHessian(hess, ridge=control$ridge)
    method_details$opg <- list(
      note='OPG-style information estimate; local prior curvature is not represented unless it appears in score variability.'
    )
  }
  
  if(uncertainty == 'fullbootstrap'){
    message('Fitting full bootstrap refits')
    fullbootstrap <- ctOptimFullBootstrapDraws(est=est, standata=standata,
      sm=sm, n=finishsamples, cores=cores, control=control,
      verbose=verbose)
    draws <- fullbootstrap$draws
    cov <- ctOptimSafeCov(stats::cov(draws), ridge=control$ridge)
    method_details$fullbootstrap <- fullbootstrap
    method_details$fullbootstrap$draws <- NULL
  }
  
  if(uncertainty %in% c('sandwich','bootstrap')){
    message('Computing score contributions')
    scoremat <- ctOptimScoreMatrix(standata=standata, sm=sm, est=est,
      cores=cores, scores=scores)
    centered <- scale(scoremat, center=TRUE, scale=FALSE)
    meat <- crossprod(centered)
    if(uncertainty == 'sandwich'){
      cov <- cov %*% meat %*% cov
      cov <- ctOptimSafeCov(cov, ridge=control$ridge)
    } else {
      draws <- ctOptimBootstrapDraws(est=est, cov=cov, scores=scoremat,
        n=finishsamples)
      cov <- ctOptimSafeCov(stats::cov(draws), ridge=control$ridge)
    }
  }
  
  if(uncertainty == 'surrogate'){
    message('Estimating local quadratic surrogate')
    surrogate <- ctOptimSurrogateHessian(est=est, lpgFunc=lpgFunc, cov=cov,
      npoints=control$surrogateNpoints, scale=control$surrogateScale,
      ridge=control$ridge,
      profile=control$surrogateProfile,
      profileTargetDrop=control$surrogateProfileTargetDrop,
      profileMaxStep=control$surrogateProfileMaxStep,
      verbose=verbose)
    hess <- surrogate$hessian
    cov <- ctOptimCovFromHessian(hess, ridge=control$ridge)
    method_details$surrogate <- surrogate
  }
  covDiagnostics <- attr(cov, 'ctOptimCovFromHessian')
  if(!is.null(covDiagnostics)) method_details$covariance <- covDiagnostics
  
  if(!is.null(hessian)) method_details$hessian <- list(
    source='exact, by forward-mode differentiation of the reverse-mode gradient')
  list(method=uncertainty, cov=cov, hessian=if(exists('hess')) hess else NULL,
    scores=scoremat, draws=draws, base_value=base[1],
    base_gradient=base_gradient, details=method_details)
}

ctOptimUpdateTransformed <- function(fit, samples, cores=1){
  savesubjectmatrices <- fit$standata$savesubjectmatrices
  sdat <- fit$standata
  if(!as.logical(savesubjectmatrices)) sdat <- standatact_specificsubjects(sdat, 1)
  fit$stanfit$transformedpars <- stan_constrainsamples(sm=fit$stanmodel,
    standata=sdat, savesubjectmatrices=savesubjectmatrices,
    savescores=fit$standata$savescores,
    dokalman=as.logical(savesubjectmatrices), samples=samples,
    cores=cores, quiet=TRUE)
  sds <- try(suppressWarnings(sqrt(diag(fit$stanfit$cov))), silent=TRUE)
  if('try-error' %in% class(sds)) sds <- rep(NA_real_, length(fit$stanfit$rawest))
  smf <- stan_reinitsf(fit$stanmodel, fit$standata)
  fit$stanfit$transformedpars_old <- NA
  try(fit$stanfit$transformedpars_old <- cbind(
    unlist(rstan::constrain_pars(smf, upars=fit$stanfit$rawest - 1.96 * sds)),
    unlist(rstan::constrain_pars(smf, upars=fit$stanfit$rawest)),
    unlist(rstan::constrain_pars(smf, upars=fit$stanfit$rawest + 1.96 * sds))),
    silent=TRUE)
  try(colnames(fit$stanfit$transformedpars_old) <- c('2.5%','mean','97.5%'),
    silent=TRUE)
  fit
}

ctOptimDataLpgFunc <- function(sm, standata, cores=1){
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  cores <- min(cores, length(unique(standata$subject)))
  
  if(cores <= 1){
    smuse <- sm
    if(!is.null(standata$recompile) && standata$recompile == 0) {
      smuse <- utils::getFromNamespace("stanmodels", "ctsem")$ctsm
    }
    smf <- stan_reinitsf(smuse, standata)
    lpg <- function(parm){
      out <- try(rstan::log_prob(smf, upars=parm, adjust_transform=TRUE,
        gradient=TRUE), silent=FALSE)
      if('try-error' %in% class(out) || is.nan(out)) {
        out <- -1e100
        attributes(out) <- list(gradient=rep(0, length(parm)))
      }
      out
    }
    return(list(lpg=lpg, cl=NULL, standata=standata, cores=1L))
  }
  
  smfile <- ''
  if(standata$recompile > 0){
    smfile <- file.path(tempdir(), paste0('ctsem_sm_',
      ceiling(stats::runif(1, 0, 100000)), '.rda'))
    save(sm, file=smfile, eval.promises=FALSE, precheck=FALSE)
  }
  cl <- makeClusterID(cores)
  parallelStanSetup(cl=cl, standata=standata, split=TRUE,
    smfile=if(standata$recompile > 0) smfile else '')
  lpg <- function(parm){
    clusterIDexport(cl, 'parm')
    out2 <- parallel::clusterEvalQ(cl=cl, parlp(parm))
    out <- try(sum(unlist(out2)), silent=TRUE)
    for(i in seq_along(out2)){
      if(i == 1) attributes(out)$gradient <- attributes(out2[[1]])$gradient
      if(i > 1) attributes(out)$gradient <-
          attributes(out)$gradient + attributes(out2[[i]])$gradient
    }
    if('try-error' %in% class(out) || is.nan(out)) {
      out <- -1e100
      attributes(out) <- list(gradient=rep(0, length(parm)))
    }
    out
  }
  list(lpg=lpg, cl=cl, standata=standata, cores=cores, smfile=smfile)
}

ctOptimFitLpgFunc <- function(fit, cores=1){
  standata <- fit$standata
  standata$savesubjectmatrices <- 0L
  ctOptimDataLpgFunc(sm=fit$stanmodel, standata=standata, cores=cores)
}

# Defined above that block, not between it and the function it documents.
# Roxygen attaches a block to whatever definition follows it, so sitting
# below it cost `ctOptimUncertainty` its export: a regenerated NAMESPACE
# carried `export(.ctResolveDraws)` in its place, un-exporting a documented
# user-facing function that three error messages tell people to call.
# How samples are produced, given how the covariance was estimated.
#
# `uncertainty` and `draws` are not independent, and pretending otherwise is
# what made this confusing: `uncertainty='is'` *is* importance sampling, so
# `draws='normal'` alongside it was accepted, warned about, and then overridden
# -- the argument appeared to be a choice and was not. It is now derived, and an
# incompatible request is an error rather than a warning about something the
# code went on to ignore.
.ctResolveDraws <- function(uncertainty, draws) {
  implied <- switch(uncertainty,
    is = 'imis',
    bootstrap = 'empirical',
    fullbootstrap = 'empirical',
    'normal')
  if(identical(draws, 'auto')) return(implied)
  if(!identical(draws, implied)) {
    stop("uncertainty='", uncertainty, "' produces draws by '", implied,
      "', so draws='", draws, "' cannot be honoured. Pass draws='auto' (the ",
      "default) or draws='", implied, "'.", call.=FALSE)
  }
  draws
}

# uncertainty='stored': redraw from the covariance the fit already carries.
#
# Every other method here builds a covariance out of model evaluations -- 2*npar
# log-probability/gradient calls for the finite-difference Hessian, more for the
# score and bootstrap methods -- and then draws from it. Asking for a different
# number of draws is not asking for that work again: the draws are iid from a
# covariance that has not changed, so this path costs no model evaluations at
# all on either backend.
#
# It replaces `ctFitAddSamples()`, which did the same thing on stan alone and
# appended rather than replaced. Appending is what makes it stan-only in spirit
# as well as in code: the rows already on the fit may have come from `is` or
# `bootstrap`, and adding normal draws to those leaves a posterior that is part
# one distribution and part another, with nothing recording the mixture.
# Replacing them cannot do that, and the warning below says when the previous
# draws were of a kind this cannot reproduce.
.ctOptimStoredRedraw <- function(fit, finishsamples, cores, verbose=0){
  julia <- inherits(fit, 'ctJuliaFit')
  cov <- if(julia) fit$estimate$cov else fit$stanfit$cov
  est <- if(julia) as.numeric(fit$estimate$raw) else fit$stanfit$rawest
  if(is.null(cov) || !length(cov) || any(!is.finite(cov))) stop(
    "uncertainty='stored' redraws from the covariance already on the fit, and ",
    "this fit has no usable one. Run ctOptimUncertainty() with a method that ",
    "computes one first -- 'hessian' is the default.", call.=FALSE)
  if(is.null(est) || length(est) != ncol(cov)) stop(
    "uncertainty='stored' needs the fit's raw estimate and stored covariance ",
    "to describe the same parameters; they do not.", call.=FALSE)

  uncertaintyfit <- if(julia) fit$uncertainty else fit$stanfit$uncertainty
  if(is.null(uncertaintyfit)) uncertaintyfit <- list(method='stored', cov=cov)
  previousdraws <- uncertaintyfit$settings$draws
  if(!is.null(previousdraws) && !identical(previousdraws, 'normal')) {
    warning("This fit's draws came from '", previousdraws, "'; redrawing from ",
      "the stored covariance gives normal draws instead. Rerun with ",
      "uncertainty='", uncertaintyfit$settings$method, "' to keep them.",
      call.=FALSE)
  }

  samples <- ctOptimNormalDraws(est, cov, finishsamples)

  # `method` keeps naming the method that produced the covariance, because that
  # is what every reader of it wants to know and it has not changed; `redrawn`
  # records that the draws were regenerated from it afterwards.
  uncertaintyfit$draws <- 'normal'
  uncertaintyfit$settings$method <- .ctJuliaOr(uncertaintyfit$settings$method,
    'stored')
  uncertaintyfit$settings$draws <- 'normal'
  uncertaintyfit$settings$finishsamples <- finishsamples
  uncertaintyfit$settings$cores <- cores
  uncertaintyfit$settings$redrawn <- TRUE

  if(julia){
    fit$estimate$rawposterior <- samples
    fit <- .ctFitNameRawUncertainty(fit)
    fit$uncertainty <- uncertaintyfit
    fit$transformedpars <- .ctBackendConstrain(fit)
    return(fit)
  }
  fit$stanfit$rawposterior <- samples
  fit <- .ctFitNameRawUncertainty(fit)
  fit$stanfit$uncertainty <- uncertaintyfit
  if(verbose > 0) message('Redrawing ', finishsamples,
    ' samples from the stored covariance')
  ctOptimUpdateTransformed(fit, samples=samples, cores=cores)
}

#' Update optimized ctsem uncertainty estimates
#'
#' Recomputes the approximate raw-parameter uncertainty for an optimized
#' \code{\link{ctFit}} object and refreshes the approximate raw-parameter
#' samples. This is the entry point for both backends; \code{ctFit} itself
#' calls it to finish an optimized fit.
#'
#' The draws it writes to \code{$rawposterior} are \emph{pseudo-posterior}
#' draws: a sample from a covariance fitted to the log-posterior surface around
#' the optimum, not a sample from the posterior itself. \code{\link{ctSample}}
#' is the other thing, genuine posterior draws by Hamiltonian Monte Carlo from
#' an optimized \code{ctJuliaFit}. Everything downstream reads either from the
#' same slot, so the difference is recorded rather than visible in the shape of
#' the result: an optimized fit carries \code{$uncertainty$settings}, a sampled
#' one carries \code{$sample}, and \code{ctOptimUncertainty} refuses a sampled
#' fit rather than replacing its draws.
#'
#' To change only the number of draws, use \code{uncertainty='stored'}, which
#' redraws from the covariance the fit already carries and evaluates no model.
#'
#' @section Backend differences:
#' The methods are the same on both backends and the covariance they produce is
#' comparable, but three things differ and are worth knowing before comparing
#' output.
#'
#' \emph{Where the result is stored.} A \code{ctStanFit} comes back with
#' \code{fit$stanfit$cov}, \code{fit$stanfit$rawposterior} and
#' \code{fit$stanfit$uncertainty}; a \code{ctJuliaFit} with
#' \code{fit$estimate$cov}, \code{fit$estimate$se},
#' \code{fit$estimate$rawposterior} and \code{fit$uncertainty}. Both record the
#' resolved settings in \code{$uncertainty$settings}. On both backends the
#' draws, the covariance and the standard errors are labelled by raw parameter,
#' so \code{fit$estimate$se['drift']} or
#' \code{apply(fit$estimate$rawposterior, 2, quantile)} reads without counting
#' columns. The two backends spell the population-SD and correlation blocks
#' differently (\code{popsd_x}/\code{rawcor_x__y} on julia,
#' \code{x_SD}/\code{y_x_corr} on stan); the order is the same.
#'
#' \emph{The Hessian.} The stan path finite-differences its gradient with a
#' single global step (\code{control$hessianStep}). The julia engine
#' differentiates its own reverse-mode gradient in forward mode, so its
#' Hessian is exact; \code{control$analyticHessian = FALSE} falls back to the
#' shared finite difference for comparison. The two agree to about 1e-4 on a
#' well-conditioned model.
#'
#' \emph{Backend-specific arguments.} \code{uncertainty='fullbootstrap'} and
#' its \code{control$bootstrapFitCores} / \code{control$bootstrapTol}, and
#' \code{control$parsteps}, are stan-only and are refused by name on a
#' \code{ctJuliaFit}. \code{control$analyticHessian} is julia-only.
#' \code{cores} means R worker processes on stan and engine threads on julia
#' (see below). The IMIS proposal defaults also differ:
#' \code{imisScaleInit = 1.1} and \code{imisTailScale = 1.1} on stan against
#' \code{1.5} and \code{1.2} on julia, which was measured -- a proposal no
#' wider than the Hessian covariance cannot correct a posterior wider than it.
#' Set them explicitly to compare the backends on this method.
#'
#' Transformed-parameter summaries are refreshed on stan and not on julia,
#' which has no parameter-matrix reconstruction through
#' \code{rstan::constrain_pars}; the raw-scale covariance and draws are
#' complete on both.
#'
#' @param fit Optimized \code{ctStanFit} or \code{ctJuliaFit} object. For a
#' \code{ctJuliaFit}, every \code{uncertainty} method except
#' \code{'fullbootstrap'} is available; that one re-optimises each resample and
#' so needs the model rebuilt rather than re-evaluated. A sampled fit of either
#' backend is refused: it already carries a posterior, and replacing it with a
#' curvature-based approximation would discard it.
#' @param uncertainty Uncertainty approximation. \code{'hessian'} uses the
#' finite-difference Hessian, \code{'surrogate'} fits a local quadratic
#' surrogate around the optimum, \code{'is'} uses Hessian-based importance
#' sampling and computes uncertainty from the weighted importance-sampling
#' distribution, \code{'bootstrap'} uses one-step score bootstrap draws with
#' Hessian bread, \code{'fullbootstrap'} resamples subjects and fully
#' re-optimizes each sample from the original maximum likelihood or MAP
#' estimate using mize L-BFGS, \code{'sandwich'} uses Hessian bread with score
#' covariance meat, and \code{'opg'} uses an OPG-style score information
#' approximation. \code{'stored'} computes nothing: it redraws
#' \code{finishsamples} normal draws from the covariance already on the fit and
#' leaves that covariance, the Hessian and the recorded \code{method}
#' untouched, marking \code{$uncertainty$settings$redrawn = TRUE}. It is the
#' way to change the number of draws, or to reseed them, without paying for the
#' covariance again -- the other methods cost at least \code{2 * npar}
#' log-probability evaluations, this one costs none -- and it warns if the
#' fit's existing draws came from \code{'is'} or \code{'bootstrap'}, which
#' normal draws from that covariance do not reproduce.
#' @param draws Approximate raw-parameter draw method. \code{'auto'} uses
#' empirical draws for \code{uncertainty='bootstrap'} and
#' \code{uncertainty='fullbootstrap'} and normal draws otherwise.
#' \code{'normal'} draws from a multivariate normal using the selected
#' covariance, \code{'empirical'} uses empirical draws when available, and
#' \code{'imis'} runs the importance sampler using the selected covariance as
#' proposal. For \code{uncertainty='is'}, \code{draws} is set to \code{'imis'}.
#' @param finishsamples Number of approximate raw-parameter samples. If
#' \code{NULL}, the existing number of rows in the fit's raw posterior
#' (\code{fit$stanfit$rawposterior} for stan, \code{fit$estimate$rawposterior}
#' for julia) is reused when available; otherwise 1000 samples are used.
#' @param cores Number of cores. If \code{NULL}, one core is used, and nothing
#' is parallelised unless a value above one is asked for. On a
#' \code{ctStanFit} these are R worker processes: each
#' log-probability/gradient evaluation is split across subjects and reassembled,
#' and score contributions and transformed quantities use them too. On a
#' \code{ctJuliaFit} there is no R cluster; the value becomes the engine's
#' subject-loop chunk ceiling for the duration of the call, and the engine caps
#' it at its own thread count. Neither route changes what is computed, though
#' both change the order things are summed in, so results are reproducible at
#' \code{cores = 1} and agree to rounding above it.
#' @param control List of method-specific options. Useful entries include
#' \code{ridge}, \code{hessianStep}, \code{surrogateNpoints},
#' \code{surrogateScale}, \code{surrogateProfile},
#' \code{surrogateProfileTargetDrop}, \code{surrogateProfileMaxStep},
#' \code{bootstrapFitCores}, \code{bootstrapTol}, \code{imisMaxIter},
#' \code{imisScaleInit}, \code{imisTailScale}, \code{isESS}, and
#' \code{isitersize}. Omitted entries use
#' \code{ridge = 1e-8}, \code{hessianStep = 1e-3},
#' \code{surrogateScale = .5}, \code{surrogateNpoints = NULL},
#' \code{surrogateProfile = TRUE},
#' \code{surrogateProfileTargetDrop = NULL},
#' \code{surrogateProfileMaxStep = 64},
#' \code{bootstrapFitCores = 1}, \code{bootstrapTol = 1e-5},
#' \code{imisMaxIter = 50}, \code{imisScaleInit = 1.1},
#' \code{imisTailScale = 1.1}, \code{isESS = 100}, and
#' \code{isitersize = 1000}. When
#' \code{surrogateNpoints} is \code{NULL}, the
#' surrogate uses at least \code{max(4 * npars, 50)} local directions. The
#' surrogate is fit in whitened coordinates relative to the proposal covariance.
#' Each direction is radially adjusted with a small evaluation budget so that
#' the retained points are close to an informative local log-probability drop,
#' rather than relying on a random cloud to land in the desired range. With
#' \code{surrogateProfile = TRUE}, all fitted surrogate curvature directions
#' are then profiled until they reach \code{surrogateProfileTargetDrop}, or the
#' surrogate target drop when \code{NULL}. The profile search uses the fitted
#' surrogate curvature to choose direction-specific starting distances, so very
#' flat directions can be checked beyond \code{surrogateProfileMaxStep} when
#' the surrogate itself predicts that a larger distance is needed. If the
#' observed profile is still flatter than the fitted surrogate predicted, a
#' small magnitude-adjusted expansion budget is used before reporting a missed
#' target. \code{parsteps} may be supplied internally to keep stepwise-fixed
#' raw parameters fixed while estimating uncertainty for the remaining
#' parameters; existing fixed indices from a previous uncertainty calculation
#' are retained when no new \code{parsteps} are supplied. It is stan-only,
#' since only \code{\link{stanoptimis}} has a stepwise phase to inherit fixed
#' parameters from, and is refused rather than ignored on a
#' \code{ctJuliaFit}. \code{analyticHessian = FALSE} is julia-only and asks for
#' the shared finite-difference Hessian instead of the engine's exact one.
#' Hessian-based covariance construction first attempts the unmodified
#' \code{solve(-hessian)} covariance and a Cholesky check. It warns when
#' numerical repair is needed, such as positive-definite projection, ridge
#' flooring of information eigenvalues, or fallback to \code{MASS::ginv};
#' diagnostics are stored in
#' \code{fit$stanfit$uncertainty$details$covariance}.
#' Score-based methods use subject-level score contributions when there are
#' at least two subjects; single-subject models warn and use case-level
#' contributions. Score-based methods warn when there are fewer than ten
#' independent subjects or no more score rows than raw parameters. Full
#' bootstrap requires at least two subjects and warns below ten independent
#' subjects. Bootstrap-style methods require at least two returned samples /
#' refits.
#' @param verbose Integer controlling progress detail.
#' @param ... Unused.
#'
#' @return The fit, of the class it came in as. The resolved method, draw
#' strategy, sample count, cores, and non-internal controls are recorded in
#' \code{fit$stanfit$uncertainty$settings} for a \code{ctStanFit} and in
#' \code{fit$uncertainty$settings} for a \code{ctJuliaFit}; see the backend
#' differences above for the other slots each writes.
#' @seealso \code{\link{ctSample}} for genuine posterior draws by Hamiltonian
#' Monte Carlo, rather than the pseudo-posterior draws from a covariance that
#' this function produces. \code{\link{ctFitAddSamples}} is the deprecated
#' stan-only predecessor of \code{uncertainty='stored'}.
#' @export
ctOptimUncertainty <- function(fit,
  uncertainty=c('hessian','surrogate','is','bootstrap','fullbootstrap',
    'sandwich','opg','stored'),
  draws=c('auto','normal','empirical','imis'), finishsamples=NULL,
  cores=NULL, control=list(), verbose=0, ...){
  
  uncertainty <- match.arg(uncertainty)
  draws <- match.arg(draws)
  # backend='julia' fits reach the same ctOptimComputeUncertainty() below,
  # through a log-probability/gradient function built from their own engine;
  # see R/ctBackendUncertainty.R.
  draws <- .ctResolveDraws(uncertainty, draws)
  if(inherits(fit, 'ctJuliaFit')) {
    # `fit$sample` (class "ctSampleDiagnostics") is set only by
    # `.ctBackendSampleAssemble()`, the routine shared by `ctFit(optimize =
    # FALSE)` and `ctSample()` -- so it marks a julia fit built from real
    # draws either way. Nothing below knows that: `.ctBackendUncertainty()`
    # treats `fit$estimate$raw` as a point estimate, builds a Hessian or score
    # matrix around it, and overwrites `fit$estimate$rawposterior` with fresh
    # curvature-based draws. Run on a sampled fit that would silently discard
    # the actual posterior draws in favour of a Gaussian approximation around
    # their mean -- wrong in a way nothing downstream would notice, since the
    # replacement is the same shape and a plausible size. The stan branch
    # below already refuses the equivalent case (`fit$stanfit$stanfit@sim`
    # populated); this mirrors it for julia.
    if(!is.null(fit$sample)) {
      stop("ctOptimUncertainty() applies to an optimized ctJuliaFit; this fit ",
        "was sampled (ctFit(optimize = FALSE) or ctSample()) and already ",
        "carries its posterior in fit$estimate$rawposterior. Read that ",
        "directly, or refit with optimize = TRUE if a curvature-based ",
        "approximation is what you want.", call.=FALSE)
    }
    # `control$parsteps` is a `stanoptimis()` concept and nothing else: the
    # stan optimiser can hold a block of raw parameters at zero for an early
    # step, and the branch below then estimates uncertainty for the remainder
    # and pads the fixed entries back in. The julia optimiser has no such
    # phase, so there is nothing for this to name -- and it reached
    # `ctOptimComputeUncertainty()`, which never reads it, as an inert list
    # element. Same standard errors as without it, and the request recorded in
    # `$uncertainty$settings$control` as though it had been honoured. Refused
    # rather than translated, because there is no julia-side meaning to
    # translate it to.
    if(!is.null(control$parsteps)) {
      stop("control$parsteps is only available for backend='stan' fits: it ",
        "holds parameters that stanoptimis() fixed during a stepwise ",
        "optimisation, and the julia optimiser has no such step. Drop it, or ",
        "refit with backend='stan'.", call.=FALSE)
    }
    # Same rule as the stan branch below, reading the julia fit's own slot.
    # Hardcoding 1000 here meant `ctOptimUncertainty(fit)` after
    # `ctOptimUncertainty(fit, finishsamples=200)` silently resampled to 1000
    # on julia and stayed at 200 on stan, against one documented default.
    if(is.null(finishsamples)) {
      finishsamples <- if(!is.null(fit$estimate$rawposterior))
        nrow(fit$estimate$rawposterior) else 1000
    }
    if(is.null(cores)) cores <- 1L
    cores <- max(1L, suppressWarnings(as.integer(cores[1])))
    if(is.na(cores)) cores <- 1L
    if(uncertainty == 'stored') return(.ctOptimStoredRedraw(fit,
      finishsamples=finishsamples, cores=cores, verbose=verbose))
    return(.ctBackendUncertainty(fit=fit, uncertainty=uncertainty, draws=draws,
      finishsamples=finishsamples, cores=cores, control=control,
      verbose=verbose))
  }
  # Named rather than asserted. "fit must be a ctStanFit object" told a caller
  # holding some other object about a class; what it needs to say is which
  # objects this does work on, since it works on both backends' fits.
  if(!'ctStanFit' %in% class(fit)) stop(
    'ctOptimUncertainty() takes an optimized ctsem fit: a ctStanFit from ',
    "ctFit(..., backend='stan') or a ctJuliaFit from ctFit(..., ",
    "backend='julia'). This object is neither.", call.=FALSE)
  if(length(fit$stanfit$stanfit@sim) > 0) {
    stop('ctOptimUncertainty currently applies to optimized ctStanFit objects')
  }
  if(is.null(finishsamples)) {
    finishsamples <- if(!is.null(fit$stanfit$rawposterior))
      nrow(fit$stanfit$rawposterior) else 1000
  }
  if(is.null(cores)) cores <- 1
  cores <- suppressWarnings(as.integer(cores[1]))
  if(!is.finite(cores) || is.na(cores) || cores < 1) cores <- 1L
  # Before `ctOptimFitLpgFunc()`, which reinitialises the stan model object,
  # and before any of it is needed: a stored redraw evaluates no model.
  if(uncertainty == 'stored') return(.ctOptimStoredRedraw(fit,
    finishsamples=finishsamples, cores=cores, verbose=verbose))
  lpg_cores <- if(uncertainty %in% c('opg','fullbootstrap') &&
      draws != 'imis') 1L else cores
  lpgsetup <- ctOptimFitLpgFunc(fit, cores=lpg_cores)
  on.exit({
    if(!is.null(lpgsetup$cl)) try(parallel::stopCluster(lpgsetup$cl), silent=TRUE)
    if(!is.null(lpgsetup$smfile) && nzchar(lpgsetup$smfile)) {
      try(file.remove(lpgsetup$smfile), silent=TRUE)
    }
  }, add=TRUE)
  if(uncertainty == 'surrogate' && is.null(control$initialCov) &&
      !is.null(fit$stanfit$cov)) {
    control$initialCov <- fit$stanfit$cov
  }
  if(is.null(control$parsteps) &&
      !is.null(fit$stanfit$uncertainty$fixedpars)) {
    control$parsteps <- fit$stanfit$uncertainty$fixedpars
  }
  parsteps <- integer()
  if(!is.null(control$parsteps)) {
    parsteps <- sort(unique(as.integer(unlist(control$parsteps))))
    parsteps <- parsteps[is.finite(parsteps) & parsteps >= 1 &
        parsteps <= length(fit$stanfit$rawest)]
    control$parsteps <- NULL
  }
  freepars <- setdiff(seq_along(fit$stanfit$rawest), parsteps)
  if(length(freepars) < 1) stop('No free parameters for uncertainty estimation')
  estuse <- fit$stanfit$rawest
  lpguse <- lpgsetup$lpg
  if(length(parsteps) > 0) {
    estuse <- fit$stanfit$rawest[freepars]
    if(!is.null(control$initialCov)) {
      control$initialCov <- control$initialCov[freepars, freepars, drop=FALSE]
    }
    lpguse <- function(parm){
      fullparm <- fit$stanfit$rawest
      fullparm[freepars] <- parm
      out <- lpgsetup$lpg(fullparm)
      grad <- attributes(out)$gradient
      if(!is.null(grad) && length(grad) == length(fullparm)) {
        attributes(out)$gradient <- grad[freepars]
      }
      out
    }
  }
  matsetup <- if(!is.null(fit$setup$matsetup)) fit$setup$matsetup else NA
  uncertaintyfit <- ctOptimComputeUncertainty(est=estuse,
    standata=lpgsetup$standata, sm=fit$stanmodel, lpgFunc=lpguse,
    uncertainty=uncertainty, finishsamples=finishsamples, cores=cores,
    matsetup=matsetup, control=control, verbose=verbose)
  if(length(parsteps) > 0) {
    freecov <- uncertaintyfit$cov
    fullcov <- diag(1e-10, length(fit$stanfit$rawest))
    fullcov[freepars, freepars] <- freecov
    uncertaintyfit$free_cov <- freecov
    uncertaintyfit$freepars <- freepars
    uncertaintyfit$fixedpars <- parsteps
    uncertaintyfit$cov <- fullcov
    if(!is.null(uncertaintyfit$hessian)) {
      fullhess <- matrix(0, length(fit$stanfit$rawest),
        length(fit$stanfit$rawest))
      fullhess[freepars, freepars] <- uncertaintyfit$hessian
      uncertaintyfit$free_hessian <- uncertaintyfit$hessian
      uncertaintyfit$hessian <- fullhess
    }
    if(!is.null(uncertaintyfit$draws)) {
      fulldraws <- matrix(rep(fit$stanfit$rawest,
          each=nrow(uncertaintyfit$draws)),
        nrow=nrow(uncertaintyfit$draws), byrow=FALSE)
      fulldraws[, freepars] <- uncertaintyfit$draws
      uncertaintyfit$free_draws <- uncertaintyfit$draws
      uncertaintyfit$draws <- fulldraws
    }
  }
  
  if(draws == 'empirical' && !is.null(uncertaintyfit$draws)) {
    samples <- uncertaintyfit$draws
  } else if(draws == 'imis'){
    if(is.null(control$imisMaxIter)) control$imisMaxIter <- 50
    if(is.null(control$imisScaleInit)) control$imisScaleInit <- 1.1
    if(is.null(control$imisTailScale)) control$imisTailScale <- 1.1
    if(is.null(control$isESS)) control$isESS <- 100
    if(is.null(control$isitersize)) control$isitersize <- 1000
    is_res <- imis_is(lpgsetup$lpg, mu_hat=fit$stanfit$rawest,
      Sigma_hat=uncertaintyfit$cov, max_iter=control$imisMaxIter,
      scale_init=control$imisScaleInit, tail_scale=control$imisTailScale,
      target_ess=control$isESS, n_batch=control$isitersize, cl=NA,
      # `verbose > 0`, not TRUE: this printed IMIS iteration progress at
      # `verbose = 0`, so the one argument meant two things across the
      # backends -- silence on julia, a page of output on stan.
      finishsamples=finishsamples, verbose=verbose > 0)
    samples <- is_res$theta
    uncertaintyfit$proposal_cov <- uncertaintyfit$cov
    if(!is.null(uncertaintyfit$details$covariance)) {
      uncertaintyfit$details$proposal_covariance <-
        uncertaintyfit$details$covariance
      uncertaintyfit$details$covariance <- NULL
    }
    if(!is.null(is_res$covariance) && all(is.finite(is_res$covariance))) {
      uncertaintyfit$cov <- ctOptimSafeCov(is_res$covariance)
    } else if(nrow(samples) > 1) {
      uncertaintyfit$cov <- ctOptimSafeCov(stats::cov(samples))
    }
    uncertaintyfit$imis <- is_res
    uncertaintyfit$details$importance_sampling <- list(ess=is_res$ess,
      df_used=is_res$df_used, covariance='weighted importance-sampling covariance')
  } else {
    samples <- ctOptimNormalDraws(fit$stanfit$rawest, uncertaintyfit$cov,
      finishsamples)
  }
  
  fit$stanfit$cov <- uncertaintyfit$cov
  fit$stanfit$rawposterior <- samples
  # As on julia. A no-op when this is called from inside `stanoptimis()`, where
  # the fit is still a stub with no model attached to read names from; `ctFit()`
  # calls it again on the assembled object.
  fit <- .ctFitNameRawUncertainty(fit)
  fit$stanfit$uncertainty <- uncertaintyfit
  fit$stanfit$uncertainty$draws <- draws
  storedControl <- control
  storedControl$initialCov <- NULL
  storedControl$parsteps <- NULL
  fit$stanfit$uncertainty$settings <- list(
    method=uncertainty,
    draws=draws,
    finishsamples=finishsamples,
    cores=cores,
    control=storedControl
  )
  if(!is.null(uncertaintyfit$scores)) fit$stanfit$subjectscores <- uncertaintyfit$scores
  message('Computing posterior approximation with ', nrow(samples), ' samples')
  fit <- ctOptimUpdateTransformed(fit, samples=samples, cores=cores)
  fit
}
