#' K fold cross validation for ctStanFit objects
#'
#' @param fit ctStanfit object
#' @param folds Number of cross validation splits to use -- 10 folds implies that the 
#' model is re-fit 10 times, each time to a data set with 1/10 of the observations randomly removed.
#' @param cores Number of processor cores to use. 
#' @param parallelFolds compute folds in parallel or use cores to finish single folds faster. 
#' parallelFolds will use folds times as much memory.
#' @param subjectwise drop random subjects instead of data rows?
#' @param keepfirstobs do not drop first observation (more stable estimates)
#' @param leaveOutN if a positive integer is given, the folds argument is ignored and 
#' instead the folds are calculated by leaving out every Nth row from the data when fitting. 
#' Leaving 2 out would result in 3 folds (starting at rows 1,2,3), each containing one third of the data. 
#' @param refit if FALSE, do not optimise parameters for the new data set, 
#' just compute the likelihoods etc from the original parameters
#' @param casewiseApproximation if TRUE, use a bootstrapped gradient contributions approach to approximate the cross validation parameters -- much faster but less reliable. 
#' @param tol tolerance for optimisation of refitted samples, can generally be more relaxed than the tolerance used for fitting initially. 
#' 
#' @details Works for \code{backend='julia'} fits as well as stan ones. There,
#' a held-out row is withheld by setting its manifest observations to missing
#' and re-preparing, which is a row the filter propagates through without an
#' update and without a likelihood contribution -- the same thing
#' \code{standata$dokalmanrows} does for the stan backend. \code{parallelFolds}
#' is ignored for julia fits, because the engine already threads its own
#' subject loop and each fold therefore uses every core.
#'
#' Be aware that each fold is an independent re-optimisation, and neither
#' backend's optimizer restarts from a flat direction. On a weakly identified
#' model, or one where withholding a fold leaves a parameter poorly determined,
#' two folds -- or two backends -- can converge to noticeably different raw
#' parameters for a similar likelihood.
#'
#' @return list
#' @export
#'
#' @examples
#' \donttest{ 
#' ctLOO(ctstantestfit)
#' }
ctLOO <- function(fit, folds = 10, cores = 2, parallelFolds = FALSE, tol = 1e-5,
  subjectwise = ifelse(length(unique(.ctFitRowSubject(fit))) >= folds, TRUE, FALSE),
  keepfirstobs = FALSE, leaveOutN = NA, refit = TRUE, casewiseApproximation = FALSE) {
  
  if(is.na(as.integer(folds))) stop('Folds must be an integer')

  # backend='julia' withholds rows at the data level rather than through
  # `standata$dokalmanrows`, and refits with the engine's own optimizer; see
  # .ctBackendLOO at the foot of this file.
  if(inherits(fit, 'ctJuliaFit')) return(.ctBackendLOO(fit=fit, folds=folds,
    cores=cores, tol=tol, subjectwise=subjectwise, keepfirstobs=keepfirstobs,
    leaveOutN=leaveOutN, refit=refit,
    casewiseApproximation=casewiseApproximation, parallelFolds=parallelFolds))

  if (!'ctStanFit' %in% class(fit) || !length(fit$stanfit$stanfit@sim) == 0) 
    stop('Not an optimized ctStanFit object')
  
  message('Using ', cores, '/', parallel::detectCores(), ' available CPU cores')
  if(all(is.na(leaveOutN))){
    if(!subjectwise) samplerows <- which(c(1, diff(fit$standata$subject)) < ifelse(keepfirstobs, 1, 999)) 
    if(subjectwise) samplerows <- 1:fit$standata$nsubjects
    
    samplerows <- sample(samplerows, length(samplerows), replace = FALSE)
    samplerows <- split(samplerows, sort(1:length(samplerows) %% folds))
    if(subjectwise) samplerows <- lapply(samplerows, function(x) which(fit$standata$subject %in% x))
  }
  
  if(!is.na(leaveOutN[1])){ 
    samplerows <- which(c(1, diff(fit$standata$subject)) < ifelse(keepfirstobs, 1, 999)) 
    samplerows <- lapply(1:(leaveOutN + 1), function(seqstart) samplerows[seq(seqstart, length(samplerows), leaveOutN + 1)])
    folds <- length(samplerows)
    message(paste0(length(samplerows), ' folds prepared...'))
  } 
  
  sdat <- fit$standata
  smodel <- fit$stanmodel
  init <- fit$stanfit$rawest
  
  # Precompute casewise gradients if subjectwise LOO is used
  if (subjectwise && is.null(fit$stanfit$subjectscores)) {
    fit$stanfit$subjectscores <- t(scorecalc(
      standata = sdat,
      est = init,
      stanmodel = smodel,
      subjectsonly = subjectwise,
      returnsubjectlist = TRUE,
      cores = cores
    ))
  }
  
  if(subjectwise) scores <- fit$stanfit$subjectscores

  
  # Parallel setup
  if (parallelFolds && cores > 1) {
    # rscript_libs: without it a worker searches its own default .libPaths(),
    # not the caller's, so library(ctsem) below can silently load a different
    # install than the one running this code (e.g. a stale globally-installed
    # package while a development tree is under test) -- the same failure
    # mode makeClusterID() (stanoptimis.R) was fixed for.
    clctsem <- parallelly::makeClusterPSOCK(min(cores, folds), rscript_libs = .libPaths())
    on.exit({ parallel::stopCluster(clctsem) }, add = TRUE)
    
    # Export required data and load library once per worker
    parallel::clusterExport(clctsem, c('sdat', 'smodel', 'init', 'parallelFolds', 'casewiseApproximation', 'scores'), envir = environment())
    parallel::clusterEvalQ(clctsem, library(ctsem))
  } else {
    clctsem <- NA  # No parallel clusters if parallelFolds is FALSE
  }
  
  # Chunk samplerows for efficiency
  chunk_size <- max(1, length(samplerows) %/% (cores * 4))  # Adjust chunk size dynamically
  samplerowssplit <- split(
    samplerows,
    rep(seq_len(ceiling(length(samplerows) / chunk_size)), each = chunk_size, length.out = length(samplerows))
  )

  # Process folds in chunks
  folded <- flexlapply(clctsem, X = samplerowssplit, fn = function(row_chunk) {
    lapply(row_chunk, function(x) {
      sdat_local <- sdat  # Avoid modifying global sdat
      if (all(is.na(leaveOutN)) || refit) {
        sdat_local$dokalmanrows[x] <- 0L
        sdat_local$dokalmanrows[-x] <- 1L
      } else {
        sdat_local$dokalmanrows[-x] <- 0L
        sdat_local$dokalmanrows[x] <- 1L
      }
      
      if(refit){ #otherwise just compute conditional on originally estimated parameters
        if (subjectwise|| casewiseApproximation) { #either use as inits or final values
          subjectsx <- unique(sdat_local$subject[x])
          score_sum <- colSums(scores[subjectsx, , drop = FALSE])
          
          # # Fisher Information approximation
          # fisher_info <- (scores) %*% t(scores) / nrow(scores)  # Empirical Fisher Information
          # fisher_info_inv <- solve(fisher_info)
          # 
          # # Parameter adjustment
          # delta_params <- -fisher_info_inv %*% score_sum

          delta_params <- -fit$stanfit$cov %*% score_sum
          init=init + delta_params
          if(casewiseApproximation) e <- list(rawest= init)
        }
        if(!casewiseApproximation) {
          e <- try(stanoptimis(standata = sdat_local,sm = smodel,init = init,
            stochastic=FALSE,carefulfit=FALSE,tol=tol,
            estonly = TRUE,cores=ifelse(parallelFolds,1,cores),plot=0,verbose=0))
        }
      }#end if refit
      if(!refit) e <- list(rawest=init) #if not refitting, just compute out of sample likelihood using estimated parameters
      if('try-error' %in% class(e)) return(NA) else{
        # sdat$savescores <- 1L
        if(all(is.na(leaveOutN)) || refit) sdat$dokalmanrows <- fit$standata$dokalmanrows #if leaveOutN and not refitting, the in sample likelihoods are what we want!
        smf <- stan_reinitsf(smodel,sdat)
        cp = rstan::constrain_pars(smf,e$rawest)
        lp = rstan::log_prob(smf,e$rawest)
        return(list(llrow=cp$llrow,logprob=lp,pars=e$rawest))
      }
    })
  }, cores = ifelse(parallelFolds, cores, 1))
  
  folded <- do.call(c, folded)
  
  llrowoos=unlist(lapply(1:folds,function(x) folded[[x]]$llrow[samplerows[[x]] ]))
  llrowoos=llrowoos[match(1:(fit$standata$ndatapoints),unlist(samplerows))]
  llrowoosSubject=sapply(unique(fit$standata$subject),function(x) 
    sum(llrowoos[fit$standata$subject==x],na.rm=FALSE) )
  llrow=fit$stanfit$transformedparsfull$llrow
  llrowSubject=sapply(unique(fit$standata$subject),function(x) 
    sum(llrow[fit$standata$subject==x],na.rm=TRUE) )
  
  llrowFolds <- as.data.table(lapply(folded, function(x) data.table(ll=c(x$llrow))[ll==0,ll:=NA]))
  
  # plot(llrow,llrowoos,col=fit$standata$subject,pch=16)
  # abline(b = 1,a=0)
  # plot(density(llrow-llrowoos),main='Original - OOS LogLik Difference')
  # abline(v=mean(llrow-llrowoos))
  
  out <- list(
    foldrows=samplerows,
    foldpars = as.matrix(data.frame(lapply(folded,function(x) x$pars))),
    # outsampleLogLikFolds=lloos,
    insampleLogLikRow=llrow,
    LogLikRowFolds = llrowFolds,
    outsampleLogLikRow=llrowoos,
    insampleLogLik=sum(llrow,na.rm=TRUE),
    outsampleLogLik=sum(llrowoos,na.rm=TRUE),
    
    insampleRowwiseEntropy = -sum(llrow,na.rm=TRUE)/fit$standata$ndatapoints,
    outsampleRowwiseEntropy = -sum(llrowoos,na.rm=TRUE)/fit$standata$ndatapoints,
    
    insampleSubjectwiseEntropy = -sum(llrow,na.rm=TRUE)/length(unique(fit$standata$subject)),
    outsampleSubjectwiseEntropy = -sum(llrowoos,na.rm=TRUE)/length(unique(fit$standata$subject)),
    
    insampleRowwiseLogLikSD = sd(llrow,na.rm=TRUE),
    outsampleRowwiseLogLikSD = sd(llrowoos,na.rm=TRUE),
    insampleSubjectwiseLogLikSD =  sd(llrowSubject,na.rm=TRUE),
    outsampleSubjectwiseLogLikSD =  sd(llrowoosSubject,na.rm=TRUE)
  )

  if(!all(is.na(leaveOutN))) out$insampleLogProb=sum(sapply(folded,function(x) x$logprob))
  
  return( out)
}






# K-fold cross validation for backend='julia' fits ---------------------------
#
# The Stan path withholds rows by zeroing `standata$dokalmanrows`, refits with
# `stanoptimis`, and reads `llrow` back out of `constrain_pars`. None of those
# three exist here, but all three have exact equivalents, and the middle one is
# the only one that needed anything new:
#
#   withhold a row  ->  set its manifests to NA and re-prepare. A fully missing
#                       row is one the filter propagates through without an
#                       update and without a likelihood contribution, which is
#                       precisely what `dokalmanrows = 0` means. It is also what
#                       `removeObs` already does for prediction, so it is a
#                       tested path rather than a new one.
#   refit           ->  `.ctJuliaOptimise()`, the same call `ctFit()` makes.
#   row likelihoods ->  the filter's own `llrow`, which the Kalman trace
#                       produces as a byproduct of the forward pass.
#
# Withholding at the *data* level rather than by a flag is what makes this
# short, and it is also stricter: there is no way for a withheld row to leak
# into the likelihood, because the engine never sees its value.
#
# What is deliberately not reproduced: `parallelFolds`. The engine threads its
# own subject loop, so folds run one at a time with every thread on each, rather
# than spawning R workers that would each need their own Julia process.

.ctBackendLOOFolds <- function(rowsubject, folds, subjectwise, keepfirstobs,
  leaveOutN) {
  nsubjects <- length(unique(rowsubject))
  # `c(1, diff(subject)) < 1` keeps only rows that are not a subject's first;
  # `< 999` keeps everything. Same expression as the Stan path, so the two
  # produce the same folds from the same seed.
  eligible <- which(c(1, diff(rowsubject)) < ifelse(keepfirstobs, 1, 999))
  if (!all(is.na(leaveOutN))) {
    out <- lapply(seq_len(leaveOutN + 1L), function(start) {
      eligible[seq(start, length(eligible), leaveOutN + 1L)]
    })
    return(out)
  }
  units <- if (subjectwise) seq_len(nsubjects) else eligible
  units <- sample(units, length(units), replace = FALSE)
  out <- split(units, sort(seq_along(units) %% folds))
  if (subjectwise) out <- lapply(out, function(x) which(rowsubject %in% x))
  out
}

# Each row's log likelihood at a given raw parameter vector, from the filter.
.ctBackendRowLoglik <- function(fit, pars) {
  spec <- .ctBackendAsModel(.ctBackendSpec(fit))
  as.numeric(.ctBackendKalmanRaw(spec, pars, subjectmatrices = FALSE,
    fields = "llrow")$llrow)
}

.ctBackendLOO <- function(fit, folds, cores, tol, subjectwise, keepfirstobs,
  leaveOutN, refit, casewiseApproximation, parallelFolds) {
  if (isTRUE(parallelFolds)) {
    warning("parallelFolds is ignored for backend='julia': the engine threads ",
      "its own subject loop, so each fold already uses every core.", call. = FALSE)
  }
  spec <- .ctBackendSpec(fit)
  model <- .ctFitModelObject(fit)
  rowsubject <- .ctFitRowSubject(fit)
  ndatapoints <- length(rowsubject)
  subjects <- unique(rowsubject)
  est <- as.numeric(fit$estimate$raw)

  message("Using ", cores, "/", parallel::detectCores(), " available CPU cores")
  samplerows <- .ctBackendLOOFolds(rowsubject, folds, subjectwise, keepfirstobs,
    leaveOutN)
  if (!all(is.na(leaveOutN))) {
    folds <- length(samplerows)
    message(folds, " folds prepared...")
  }

  # Score-based approximation, when asked for: one Newton step away from the
  # estimate using the held-out subjects' own score contributions and the
  # fitted covariance. Same construction as the Stan path, and the engine
  # produces the scores directly rather than by re-initialising per subject.
  scores <- NULL
  if (isTRUE(casewiseApproximation) || (subjectwise && isTRUE(refit))) {
    scores <- try(.ctBackendScoreMatrix(fit, est), silent = TRUE)
    if (inherits(scores, "try-error")) scores <- NULL
  }
  covariance <- fit$estimate$cov

  folded <- lapply(seq_len(folds), function(foldi) {
    heldout <- samplerows[[foldi]]
    pars <- est
    if (isTRUE(refit)) {
      start <- est
      if (!is.null(scores) && !is.null(covariance)) {
        held <- unique(rowsubject[heldout])
        held <- held[held <= nrow(scores)]
        if (length(held)) {
          start <- as.numeric(est - covariance %*% colSums(scores[held, , drop = FALSE]))
        }
      }
      if (isTRUE(casewiseApproximation)) {
        pars <- start
      } else {
        training <- spec$data
        # `[-heldout, ]` would drop the rows; setting them missing keeps the
        # subject's timeline intact, so the filter still propagates across the
        # gap exactly as it will when the row is scored out of sample.
        training[heldout, model$manifestNames] <- NA
        trainingfit <- .ctFitReplaceData(fit, training)
        result <- try(.ctJuliaOptimise(trainingfit$model_spec, start,
          backendcontrol = fit$args$resolved$backendcontrol, cores = cores, tol = tol),
          silent = TRUE)
        if (inherits(result, "try-error")) return(NULL)
        pars <- as.numeric(result$minimizer)
      }
    }
    # Scored against the *full* data, so the held-out rows have a likelihood to
    # report; the Stan path restores `dokalmanrows` before this step for the
    # same reason. When `leaveOutN` is used without refitting, the in-sample
    # likelihoods are the ones wanted, and they are the same call.
    list(llrow = .ctBackendRowLoglik(fit, pars), pars = pars)
  })

  usable <- !vapply(folded, is.null, logical(1L))
  if (!any(usable)) stop("Every cross-validation fold failed to refit.", call. = FALSE)
  if (any(!usable)) {
    warning(sum(!usable), " of ", folds, " folds failed to refit and are dropped.",
      call. = FALSE)
  }

  llrowoos <- rep(NA_real_, ndatapoints)
  for (foldi in which(usable)) {
    rows <- samplerows[[foldi]]
    llrowoos[rows] <- folded[[foldi]]$llrow[rows]
  }
  llrowoosSubject <- vapply(subjects, function(s)
    sum(llrowoos[rowsubject == s], na.rm = FALSE), numeric(1L))

  llrow <- .ctBackendRowLoglik(fit, est)
  llrow[llrow == 0] <- NA
  llrowSubject <- vapply(subjects, function(s)
    sum(llrow[rowsubject == s], na.rm = TRUE), numeric(1L))

  llrowFolds <- data.table::as.data.table(lapply(folded[usable], function(x) {
    value <- x$llrow
    value[value == 0] <- NA
    value
  }))

  out <- list(
    foldrows = samplerows,
    foldpars = do.call(cbind, lapply(folded[usable], function(x) x$pars)),
    # A 1-row matrix, not a vector: that is the shape the Stan path returns
    # (it comes straight out of `transformedparsfull$llrow`), and code that
    # subsets it as `[1, ]` should not care which backend produced the fit.
    insampleLogLikRow = matrix(llrow, nrow = 1L),
    LogLikRowFolds = llrowFolds,
    outsampleLogLikRow = llrowoos,
    insampleLogLik = sum(llrow, na.rm = TRUE),
    outsampleLogLik = sum(llrowoos, na.rm = TRUE),

    insampleRowwiseEntropy = -sum(llrow, na.rm = TRUE) / ndatapoints,
    outsampleRowwiseEntropy = -sum(llrowoos, na.rm = TRUE) / ndatapoints,

    insampleSubjectwiseEntropy = -sum(llrow, na.rm = TRUE) / length(subjects),
    outsampleSubjectwiseEntropy = -sum(llrowoos, na.rm = TRUE) / length(subjects),

    insampleRowwiseLogLikSD = stats::sd(llrow, na.rm = TRUE),
    outsampleRowwiseLogLikSD = stats::sd(llrowoos, na.rm = TRUE),
    insampleSubjectwiseLogLikSD = stats::sd(llrowSubject, na.rm = TRUE),
    outsampleSubjectwiseLogLikSD = stats::sd(llrowoosSubject, na.rm = TRUE)
  )
  out
}
