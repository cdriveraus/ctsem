#' Cross validation for ctsem fits: K fold by refitting, or PSIS for Laplace fits
#'
#' Estimates out-of-sample predictive performance, by refitting the model with
#' each fold of the data held out in turn or, for a julia Laplace fit, by
#' Pareto-smoothed importance sampling without refitting.
#'
#' @param fit A fit from \code{\link{ctFit}}, on either backend.
#' @param folds Number of cross validation splits to use -- 10 folds implies that the
#' model is re-fit 10 times, each time to a data set with 1/10 of the observations randomly removed.
#' @param cores Number of processor cores to use.
#' @param parallelFolds compute folds in parallel or use cores to finish single folds faster.
#' parallelFolds will use folds times as much memory.
#' @param subjectwise drop random subjects instead of data rows? With
#'   \code{method = 'psis'} this chooses the level instead -- leave one subject
#'   (unit) out when \code{TRUE}, the default there, and leave one row out when
#'   \code{FALSE}.
#' @param keepfirstobs do not drop first observation (more stable estimates)
#' @param leaveOutN if a positive integer is given, the folds argument is ignored and
#' instead the folds are calculated by leaving out every Nth row from the data when fitting.
#' Leaving 2 out would result in 3 folds (starting at rows 1,2,3), each containing one third of the data.
#' @param refit if FALSE, do not optimise parameters for the new data set,
#' just compute the likelihoods etc from the original parameters
#' @param casewiseApproximation \code{TRUE} replaces each fold's refit by one
#'   Newton step from the full estimate: the held-out data's score times the
#'   fit's covariance. Much faster than refitting, and close to it when no fold
#'   moves the estimate far. On a julia fit each fold's influence, the step's
#'   squared length in the covariance's metric, is returned in
#'   \code{foldinfluence}, and a number instead of \code{TRUE} steps only the
#'   folds whose influence is below it and refits the rest (see Details). Needs an optimized fit with
#'   a covariance; on a stan fit, subject folds and \code{TRUE} or
#'   \code{FALSE} only.
#' @param tol objective tolerance for each fold's refit. On a stan fit it
#'   defaults to 1e-5. On a julia fit a fold by default stops once a step is
#'   predicted to gain less than 1e-4 in log likelihood, which keeps each
#'   held-out subject's error to about 0.014 times the square root of its
#'   influence (see \code{casewiseApproximation}); an objective tolerance
#'   stops it while the gradient can still be large, and at 1e-5 moved one
#'   40-subject model's out-of-sample log likelihood by 1 in total.
#' @param method \code{'kfold'} (the default) refits the model once per fold.
#'   \code{'psis'} is refit-free Pareto-smoothed importance sampling
#'   leave-one-out for \code{backend = 'julia'} fits (see Details); it needs
#'   the \pkg{loo} package. The K-fold arguments \code{folds}, \code{parallelFolds},
#'   \code{tol}, \code{keepfirstobs}, \code{leaveOutN}, \code{refit} and
#'   \code{casewiseApproximation} mean nothing to it and are refused if given.
#' @param nsamples Number of importance sampling draws for \code{method =
#'   'psis'}. A sampled fit's posterior draws are used instead: all of them,
#'   or thinned to \code{nsamples} when it is given. Fewer draws estimate the
#'   Pareto k of a heavy-tailed correction too low: on one fit with an
#'   outlying subject the proposal's k read 0.54 at 1000 draws and 0.83 at
#'   4000.
#'
#' @details Works for \code{backend='julia'} fits as well as stan ones. There,
#' a held-out row is withheld by setting its manifest observations to missing
#' and re-preparing, which is a row the filter propagates through without an
#' update and without a likelihood contribution -- the same thing
#' \code{standata$dokalmanrows} does for the stan backend. For julia fits
#' \code{parallelFolds = TRUE} runs up to \code{cores} folds at once, each in a
#' background R session with a one-thread Julia, started and compiled for the
#' model before the first fold (it needs the \pkg{future} package; the
#' sessions persist for the next call, and \code{\link{ctJuliaKill}('workers')}
#' releases them). Without it the folds run one at a time on \code{cores}
#' threads. Folds are separate fits, so they gain more from processes than
#' from threads, except where each fold takes only seconds and starting the
#' sessions costs more than it saves. \code{$processes} on the result says
#' which ran.
#'
#' Be aware that each fold is an independent re-optimisation, and neither
#' backend's optimizer restarts from a flat direction. On a weakly identified
#' model, or one where withholding a fold leaves a parameter poorly determined,
#' two folds -- or two backends -- can converge to noticeably different raw
#' parameters for a similar likelihood.
#'
#' \strong{K fold on an \code{intoverpop = 'laplace'} fit.} The random effects
#' are integrated out, so what is scored depends on the fold type. With row
#' folds, a held-out row is scored with each subject's random effects at their
#' mode given the fold's \emph{training} rows, at the fold's parameters -- the
#' held-out row never informs the mode it is scored at. With subject folds a
#' held-out subject's rows are not scored at all: the out-of-sample quantity is
#' the subject's Laplace marginal likelihood at the fold's parameters,
#' returned as \code{outsampleLogLikSubject}, and \code{outsampleLogLikRow} is
#' \code{NA} because the marginal does not decompose into rows. In a multilevel
#' model a subject's marginal is conditional on the rest of its outermost group,
#' as the difference of that group's marginal with and without the held-out
#' subjects; several subjects of one group held out together share that joint
#' value evenly. \code{insampleLogLik} is correspondingly the Laplace marginal at
#' the estimate for subject folds, and the sum of the rows at the full-data
#' modes for row folds.
#'
#' \strong{\code{casewiseApproximation}.} Each fold's refit is replaced by one
#' Newton step from the full estimate, \code{est - cov \%*\% s}, where
#' \code{s} is the score of the held-out data at the estimate: each unit the
#' fold empties (a subject, or under the Laplace layer an outermost group)
#' contributes its likelihood score, and a unit it leaves partly observed -- a
#' row fold, or a subject inside a group -- the difference between its score on
#' the full and on the training data, at one gradient of the training objective
#' per fold. The prior stays, as it does in a refit. This is the one-step
#' (infinitesimal jackknife) approximation to the refit; its error grows faster
#' than the step. A fold's influence \code{t(s) \%*\% cov \%*\% s} measures the step in
#' posterior standard deviations squared, and a number given as
#' \code{casewiseApproximation} refits exactly every fold whose influence is at
#' or above it, so the approximation is used only where the step is short.
#' \code{foldstepped} says which folds were stepped. On a julia fit a stepped
#' fold scores its held-out rows only, together with every other stepped fold
#' in one pass over the data, so its column of \code{LogLikRowFolds} is
#' \code{NA} elsewhere: the training rows at each step would cost a pass a
#' fold, which for leave-one-subject-out is the whole cost.
#'
#' \strong{\code{method = 'psis'}.} Pareto-smoothed importance sampling
#' (Vehtari, Gelman and Gabry, 2017), without refitting. It estimates the
#' Bayesian leave-one-out predictive density -- a unit's likelihood averaged
#' over the posterior without it -- where K fold scores each fold at its refit
#' estimate, so the two differ by more than Monte Carlo error, most for a unit
#' that moves the estimate. Two levels:
#' \itemize{
#'   \item \code{subjectwise = TRUE}: leave one unit out, where a unit is the
#'   outermost level of random effect -- a subject in a one-level model, a
#'   group in a multilevel one, since subjects sharing a group effect are not
#'   independent given the population parameters. A unit's term is its
#'   likelihood with its random effects integrated out: by the filter under
#'   \code{intoverpop = 'augmented'}, by the Laplace approximation otherwise.
#'   On a sampled fit (\code{optimize = FALSE}, or
#'   \code{ctFitUncertainty(fit, 'sample')}) the population parameters are
#'   its posterior draws. On an optimized fit they are drawn from its normal
#'   approximation (\code{fit$estimate$cov}, so the fit must carry one),
#'   widened by a factor of 1.5 and importance corrected to the posterior. A
#'   unit whose Pareto k is above 0.7 is estimated again with \code{nsamples}
#'   more draws, centred where leaving it out moves the estimate -- the step
#'   \code{casewiseApproximation} takes -- pooled with draws around the
#'   estimate, and the estimate with the lower k is kept
#'   (\code{pointwise$refined}). That needs the fit's objective to be the
#'   posterior it describes, so it is not done for a fit sampled jointly over
#'   parameters and random effects (\code{intoverpop = 'none'}), whose unit
#'   terms are Laplace approximations to the sampled marginal. When the normal
#'   approximation itself misses the posterior (\code{proposal_k}, printed,
#'   above 0.7) every unit's estimate is unreliable, and a warning says to
#'   sample the fit and use its draws instead. No warning is not proof the
#'   approximation is good: k is estimated too low from few draws, so a value
#'   approaching 0.7, or a posterior you expect to be skewed (an outlying
#'   subject, a variance near zero), is worth rerunning with more
#'   \code{nsamples} or on a sampled fit.
#'   \item \code{subjectwise = FALSE}, on an optimized
#'   \code{intoverpop = 'laplace'} fit only: leave one row out, conditional on the
#'   estimated population parameters and integrating over each unit's random
#'   effects, drawn from the Gaussian the Laplace approximation fits at their
#'   mode, widened by a factor of 1.5. This is pointwise leave-one-out of the
#'   filter's one-step-ahead terms -- each row's density given the rows before
#'   it -- not exact time series
#'   leave-one-out, and it does not propagate uncertainty in the population
#'   parameters.
#' }
#' Draws whose inner mode solve fails are dropped and counted. A Pareto k above
#' 0.7 flags a unit or row whose estimate is unreliable; K fold is then the
#' better tool for it.
#'
#' \strong{Both methods cross-validate the Laplace objective.} A Laplace fit's
#' estimate is corrected by quadrature by default
#' (\code{optimcontrol$laplace_correct}, see \code{\link{ctFit}}). The folds are
#' refits of the Laplace objective, scored by it, and are not corrected -- each
#' would need its own Hessian and quadrature gradient -- so the in-sample side,
#' the folds' starting point and the PSIS proposal all use the Laplace optimum,
#' \code{fit$laplace$correction$laplace_estimate}, rather than the corrected
#' estimate. \code{scoring} says so when the fit was corrected.
#'
#' \strong{A fold's estimate is never certified.} Whatever the route, a refit
#' fold is optimised with \code{optimcontrol$certify = FALSE}: certifying and
#' correcting a fold the way the full fit is would cost each fold its own
#' Hessian, which is the expense K-fold exists to avoid. \code{scoring} says so
#' for any refit (\code{refit = TRUE}, the default, and not
#' \code{casewiseApproximation}), on both routes.
#'
#' @return For \code{method = 'kfold'}, a list with \code{foldrows},
#'   \code{foldpars}, \code{insampleLogLikRow}, \code{LogLikRowFolds},
#'   \code{outsampleLogLikRow}, \code{insampleLogLik}, \code{outsampleLogLik}
#'   and entropy and standard deviation summaries, and \code{processes},
#'   whether the folds ran in worker processes. On a julia fit with
#'   \code{casewiseApproximation}, \code{foldinfluence} and
#'   \code{foldstepped} give each fold's influence and whether it was stepped
#'   rather than refit. \code{scoring} is added
#'   whenever the folds were refit, naming what was scored and that they are
#'   uncertified (see Details); an \code{intoverpop = 'laplace'} fit also adds
#'   \code{insampleLogLikSubject} and \code{outsampleLogLikSubject}. For
#'   \code{method = 'psis'}, an object of class
#'   \code{ctLOOpsis}: \code{elpd_loo}, \code{se_elpd_loo}, \code{p_loo},
#'   \code{looic}, \code{pointwise} (a data frame of \code{elpd_loo},
#'   \code{lpd}, \code{p_loo} and \code{pareto_k} per unit or row, and at
#'   the unit level \code{refined}),
#'   \code{n_high_k}, \code{ndraws}, \code{ndropped}, \code{proposal_k} (the
#'   Pareto k of the correction from proposal to posterior: one value at the
#'   unit level, \code{NA} for a sampled fit, one per unit at the row level),
#'   \code{level} and \code{note}.
#' @references Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian
#'   model evaluation using leave-one-out cross-validation and WAIC.
#'   \emph{Statistics and Computing}, 27, 1413--1432.
#' @export
#'
#' @examples
#' \donttest{
#' ctLOO(ctstantestfit)
#' }
ctLOO <- function(fit, folds = 10, cores = getOption("mc.cores", 2L), parallelFolds = FALSE, tol = NULL,
  subjectwise = ifelse(length(unique(.ctFitRowSubject(fit))) >= folds, TRUE, FALSE),
  keepfirstobs = FALSE, leaveOutN = NA, refit = TRUE, casewiseApproximation = FALSE,
  method = c("kfold", "psis"), nsamples = 1000) {

  method <- match.arg(method)
  if (identical(method, "psis")) {
    if (!.ctFitIsJulia(fit)) {
      stop("method = 'psis' needs a backend = 'julia' fit; use method = ",
        "'kfold' for this fit.", call. = FALSE)
    }
    kfoldonly <- c("folds", "parallelFolds", "tol", "keepfirstobs", "leaveOutN",
      "refit", "casewiseApproximation")
    given <- intersect(names(match.call())[-1L], kfoldonly)
    if (length(given)) {
      stop("method = 'psis' does not refit, so ", paste(given, collapse = ", "),
        if (length(given) > 1L) " are" else " is",
        " not used; drop ", if (length(given) > 1L) "them" else "it",
        " or use method = 'kfold'.", call. = FALSE)
    }
    if (missing(subjectwise)) subjectwise <- TRUE
    return(.ctBackendLOOPsis(fit, subjectwise = subjectwise, ndraws = nsamples,
      cores = cores, thin = !missing(nsamples)))
  }
  if (!missing(nsamples)) {
    stop("nsamples is used only by method = 'psis'.", call. = FALSE)
  }

  if(is.na(as.integer(folds))) stop('Folds must be an integer')
  casewiseBar <- .ctLOOCasewiseBar(casewiseApproximation)
  # The Stan path steps by whole subjects' scores (`scorecalc()`), so a fold
  # holding out only some of a subject's rows would be approximated as one
  # holding out all of them; it stopped there on an undefined `scores`.
  if (!inherits(fit, 'ctJuliaFit') && casewiseBar > 0) {
    if (!isTRUE(casewiseApproximation)) {
      stop("A numeric casewiseApproximation, an influence bar, needs a ",
        "backend = 'julia' fit; use TRUE or FALSE here.", call. = FALSE)
    }
    if (!isTRUE(subjectwise) || !all(is.na(leaveOutN))) {
      stop("On a stan fit casewiseApproximation removes whole subjects, so it ",
        "needs subject folds: subjectwise = TRUE and no leaveOutN.", call. = FALSE)
    }
  }

  # backend='julia' withholds rows at the data level rather than through
  # `standata$dokalmanrows`, and refits with the engine's own optimizer; see
  # .ctBackendLOO at the foot of this file.
  if(inherits(fit, 'ctJuliaFit')) return(.ctBackendLOO(fit=fit, folds=folds,
    cores=cores, tol=tol, subjectwise=subjectwise, keepfirstobs=keepfirstobs,
    leaveOutN=leaveOutN, refit=refit,
    casewiseApproximation=casewiseApproximation, parallelFolds=parallelFolds))

  if (!'ctStanFit' %in% class(fit) || !length(fit$stanfit$stanfit@sim) == 0) 
    stop('Not an optimized ctStanFit object')
  if (is.null(tol)) tol <- 1e-5
  
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
    # mode makeClusterID() (stanoptimis.R) was fixed for. It narrows the gap
    # rather than closing it -- see .ctClusterCheckBuild(), called below.
    clctsem <- parallelly::makeClusterPSOCK(min(cores, folds), rscript_libs = .libPaths())
    on.exit({ parallel::stopCluster(clctsem) }, add = TRUE)
    
    # Export required data and load library once per worker
    parallel::clusterExport(clctsem, c('sdat', 'smodel', 'init', 'parallelFolds', 'casewiseApproximation', 'scores'), envir = environment())
    parallel::clusterEvalQ(clctsem, library(ctsem))
    .ctClusterCheckBuild(clctsem)
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
    # Whether the folds ran in a cluster of worker processes, as the julia
    # path reports.
    processes = isTRUE(parallelFolds && cores > 1),
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

# Whether a fit integrates its random effects by Laplace. A fit that samples
# them (`intoverpop = 'none'`) carries the same `spec$laplace` description, so
# the structure's presence is not enough; see `.ctBackendIntOverPop`.
.ctBackendIsLaplace <- function(fit) {
  if (!.ctFitIsJulia(fit)) return(FALSE)
  spec <- .ctBackendSpec(fit)
  !is.null(spec$laplace) && identical(.ctBackendIntOverPop(spec), "laplace")
}

# Each row's log likelihood at a given raw parameter vector, from the filter.
#
# `modesFrom`, for a Laplace fit, is the specification whose data the random
# effect modes are solved from: a fold's training data. Without it the modes
# come from the fit's own data, held-out rows included, and a held-out row is
# scored at a mode that has already seen it. Measured: shifting one held-out
# observation by five units moved the score of an *earlier* row of the same
# subject from -0.82 to -1.54, through the mode alone. The mechanism is the one
# `ctKalman(removeObs=)` uses; see `.ctBackendKalmanSpec`.
.ctBackendRowLoglik <- function(fit, pars, modesFrom = NULL) {
  spec <- .ctBackendAsModel(.ctBackendSpec(fit))
  if (!is.null(modesFrom)) {
    idname <- .ctFitModelObject(fit)$subjectIDname
    map <- match(unique(spec$data[[idname]]), unique(modesFrom$data[[idname]]))
    if (anyNA(map)) {
      stop("A fold's training data lost a subject, so its random effects cannot ",
        "be taken from that fold.", call. = FALSE)
    }
    attr(spec, "laplaceSource") <- .ctBackendAsModel(modesFrom)
    attr(spec, "laplaceSubjects") <- map
  }
  as.numeric(.ctBackendKalmanRaw(spec, pars, subjectmatrices = FALSE,
    fields = "llrow")$llrow)
}

# Each unit's log likelihood (`unit_loglik`, units x draws) and the log
# posterior (`value`) at each column of `draws` (`npar x ndraws`), whether each
# draw could be evaluated, and the unit each subject is in. A unit is the
# outermost level of random effect: a subject on the filter's route, which has
# no other, and a subject or group under the Laplace layer. The engine's
# `ctsem_unit_terms` has a method for each objective, so this does not ask
# which it holds. One bridge call for all columns. `subjects`, Laplace layer
# only, scores just the units holding them: no value, every other unit NA.
.ctBackendUnitTerms <- function(spec, draws, subjects = NULL) {
  draws <- as.matrix(draws)
  storage.mode(draws) <- "double"
  module <- .ctJuliaModule(spec$project)
  arguments <- list(.ctJuliaObjective(.ctBackendAsModel(spec), hybrid = TRUE),
    .ctJuliaPut(draws))
  if (length(subjects)) arguments$subjects <- .ctJuliaVector(as.integer(subjects))
  res <- .ctBackendJuliaValue(do.call(module$ctsem_unit_terms, arguments))
  list(value = as.numeric(res$value),
    unit_loglik = matrix(as.numeric(res$unit_loglik), ncol = ncol(draws)),
    converged = as.logical(res$converged), unit = as.integer(res$unit))
}

# A held-out subject's out-of-sample value in a subject fold: its unit's
# marginal on the full data less the same unit's marginal on the training data,
# both at the fold's parameters. For a one-level model the training term of a
# subject with every row withheld is exactly zero, so this is the subject's own
# marginal; nested in a group, it is the subject's marginal *given its group
# mates*, which is the quantity that leaves it out. Several held-out subjects in
# one unit get their joint value, shared evenly -- the joint does not split.
.ctBackendLOOHeldSubjects <- function(full, training, held) {
  if (!identical(full$unit, training$unit)) {
    stop("A fold's training data changed the grouping of subjects into units.",
      call. = FALSE)
  }
  out <- rep(NA_real_, length(full$unit))
  for (U in unique(full$unit[held])) {
    members <- held[full$unit[held] == U]
    out[members] <- (full$unit_loglik[U, 1L] - training$unit_loglik[U, 1L]) /
      length(members)
  }
  out
}

.ctBackendLOO <- function(fit, folds, cores, tol, subjectwise, keepfirstobs,
  leaveOutN, refit, casewiseApproximation, parallelFolds) {
  spec <- .ctBackendSpec(fit)
  model <- .ctFitModelObject(fit)
  rowsubject <- .ctFitRowSubject(fit)
  ndatapoints <- length(rowsubject)
  subjects <- unique(rowsubject)
  # The Laplace optimum, not a quadrature-corrected estimate. Every fold is a
  # refit of the Laplace objective and is scored by it, so the in-sample side
  # has to be that objective's optimum too: a corrected full-data estimate
  # against uncorrected folds would compare two different estimators. The
  # folds are not corrected themselves -- each would need a Hessian and
  # `2 * npar` quadrature evaluations -- so what is cross-validated is the
  # Laplace-approximated model, and `scoring` says so.
  est <- .ctLaplaceOptimum(fit)
  corrected <- .ctLaplaceIsCorrected(fit)

  message("Using ", cores, "/", parallel::detectCores(), " available CPU cores")
  samplerows <- .ctBackendLOOFolds(rowsubject, folds, subjectwise, keepfirstobs,
    leaveOutN)
  if (!all(is.na(leaveOutN))) {
    folds <- length(samplerows)
    message(folds, " folds prepared...")
  }

  # A Laplace fit's random effects are integrated out, so the fold's training
  # rows have to be kept away from them as well as from the parameters. The
  # condition matches the stan path's: the rows are withheld unless `leaveOutN`
  # is used without refitting, where the in-sample likelihoods are the ones
  # wanted. Subject folds on such a fit score each held-out subject's marginal
  # rather than its rows; see the Details of ?ctLOO.
  laplace <- .ctBackendIsLaplace(fit)
  withheld <- all(is.na(leaveOutN)) || isTRUE(refit)
  subjectlevel <- laplace && isTRUE(subjectwise) && all(is.na(leaveOutN))

  # In sample, at the estimate. A fold that is not refit and takes no modes
  # from its training rows scores exactly this, so it reads it from here.
  insample <- .ctBackendRowLoglik(fit, est)
  # Each unit's term at the estimate, and the unit each subject is in.
  fitted <- if (subjectlevel) .ctBackendUnitTerms(spec, est) else NULL

  # One Newton step from the estimate by the held-out data's score
  # (`.ctLOOFoldScore()`): a fold's answer under casewiseApproximation, below
  # its influence bar, and otherwise a refit's start. The score rows are the
  # units -- subjects on the filter's route, outermost groups under the
  # Laplace layer -- and carry no prior, which a refit keeps.
  casewiseBar <- .ctLOOCasewiseBar(casewiseApproximation)
  stepping <- casewiseBar > 0 && isTRUE(refit)
  if (stepping && !is.null(fit$sample)) {
    stop("casewiseApproximation approximates refitting an optimized fit, and ",
      "this fit was sampled; use method = 'psis', which uses its posterior ",
      "draws.", call. = FALSE)
  }
  scores <- NULL
  unitOf <- NULL
  covariance <- fit$estimate$cov
  if (isTRUE(refit)) {
    scores <- try(.ctBackendScoreMatrix(fit, est, prior = FALSE), silent = TRUE)
    if (stepping && (inherits(scores, "try-error") || is.null(covariance))) {
      stop("casewiseApproximation needs the held-out data's scores and the ",
        "fit's covariance, and ", if (is.null(covariance)) {
          "this fit has no covariance (fit with estonly = FALSE, or run ctFitUncertainty())"
        } else conditionMessage(attr(scores, "condition")),
        ". Use casewiseApproximation = FALSE.", call. = FALSE)
    }
    if (inherits(scores, "try-error")) {
      scores <- NULL
    } else {
      unitOf <- if (!is.null(fitted)) fitted$unit else if (!is.null(spec$laplace)) {
        .ctBackendUnitTerms(spec, est, subjects = 1L)$unit
      } else seq_along(subjects)
    }
  }

  # Everything a fold reads, in one list, so a fold can run in a worker
  # process as well as here (`.ctBackendLOOFold()`).
  context <- list(fit = fit, spec = spec, model = model, est = est,
    rowsubject = rowsubject, refit = refit, laplace = laplace,
    withheld = withheld, subjectlevel = subjectlevel, scores = scores,
    unitOf = unitOf, covariance = covariance, stepping = stepping,
    casewiseBar = casewiseBar, tol = tol, cores = cores, insample = insample)
  jobs <- function(context) lapply(seq_len(folds), function(foldi)
    list(heldout = samplerows[[foldi]], context = context))
  # `parallelFolds`: each fold whole in a worker process with a one-thread
  # Julia (`.ctBackendWorkerMap()`), up to `cores` at a time, rather than one
  # fold at a time on every thread. Folds are separate fits, and separate fits
  # spread over processes rather than threads (see the worker map). Falls back
  # to this session when the pool cannot be used.
  folded <- NULL
  if (isTRUE(parallelFolds) && folds > 1L) {
    context$cores <- 1L
    folded <- .ctBackendWorkerMap(fit, jobs(context), ".ctBackendLOOFold",
      workers = cores, values = est)
    context$cores <- cores
    if (is.null(folded)) message("Folds run one at a time in this session instead.")
  }
  processes <- !is.null(folded)
  if (!processes) folded <- lapply(jobs(context), .ctBackendLOOFold)

  usable <- !vapply(folded, is.null, logical(1L))
  if (!any(usable)) stop("Every cross-validation fold failed to refit.", call. = FALSE)
  if (any(!usable)) {
    warning(sum(!usable), " of ", folds, " folds failed to refit and are dropped.",
      call. = FALSE)
  }
  deferred <- which(usable)[vapply(folded[usable], function(x) isTRUE(x$deferred),
    logical(1L))]
  if (length(deferred)) {
    folded <- .ctBackendLOOScoreSteps(fit, est, folded, deferred, samplerows,
      rowsubject)
  }

  llrow <- insample
  llrow[llrow == 0] <- NA
  llrowSubject <- vapply(subjects, function(s)
    sum(llrow[rowsubject == s], na.rm = TRUE), numeric(1L))

  llrowoos <- rep(NA_real_, ndatapoints)
  if (subjectlevel) {
    # The marginal does not decompose into rows, so there is no per-row
    # out-of-sample number that means what the subject value means; the rows
    # stay NA rather than holding something else under the same name.
    llrowoosSubject <- rep(NA_real_, length(subjects))
    for (foldi in which(usable)) {
      value <- folded[[foldi]]$subject
      llrowoosSubject[!is.na(value)] <- value[!is.na(value)]
    }
    unconverged <- sum(!vapply(folded[usable], function(x) isTRUE(x$converged),
      logical(1L)))
    if (unconverged) {
      warning("The inner mode solve did not converge in ", unconverged, " of ",
        sum(usable), " folds; their subject values are approximate.", call. = FALSE)
    }
    size <- tabulate(fitted$unit, nbins = nrow(fitted$unit_loglik))
    insampleSubject <- fitted$unit_loglik[fitted$unit, 1L] / size[fitted$unit]
    insampleTotal <- sum(fitted$unit_loglik[, 1L])
    outsampleTotal <- sum(llrowoosSubject)
  } else {
    for (foldi in which(usable)) {
      rows <- samplerows[[foldi]]
      llrowoos[rows] <- folded[[foldi]]$llrow[rows]
    }
    llrowoosSubject <- vapply(subjects, function(s)
      sum(llrowoos[rowsubject == s], na.rm = FALSE), numeric(1L))
    insampleSubject <- llrowSubject
    insampleTotal <- sum(llrow, na.rm = TRUE)
    outsampleTotal <- sum(llrowoos, na.rm = TRUE)
  }

  llrowFolds <- data.table::as.data.table(lapply(folded[usable], function(x) {
    value <- x$llrow
    value[value == 0] <- NA
    value
  }))

  out <- list(
    foldrows = samplerows,
    # Whether the folds ran in worker processes (`parallelFolds`) or here.
    processes = processes,
    foldpars = do.call(cbind, lapply(folded[usable], function(x) x$pars)),
    # A 1-row matrix, not a vector: that is the shape the Stan path returns
    # (it comes straight out of `transformedparsfull$llrow`), and code that
    # subsets it as `[1, ]` should not care which backend produced the fit.
    insampleLogLikRow = matrix(llrow, nrow = 1L),
    LogLikRowFolds = llrowFolds,
    outsampleLogLikRow = llrowoos,
    insampleLogLik = insampleTotal,
    outsampleLogLik = outsampleTotal,

    insampleRowwiseEntropy = -insampleTotal / ndatapoints,
    outsampleRowwiseEntropy = -outsampleTotal / ndatapoints,

    insampleSubjectwiseEntropy = -insampleTotal / length(subjects),
    outsampleSubjectwiseEntropy = -outsampleTotal / length(subjects),

    insampleRowwiseLogLikSD = stats::sd(llrow, na.rm = TRUE),
    outsampleRowwiseLogLikSD = stats::sd(llrowoos, na.rm = TRUE),
    insampleSubjectwiseLogLikSD = stats::sd(insampleSubject, na.rm = TRUE),
    outsampleSubjectwiseLogLikSD = stats::sd(llrowoosSubject, na.rm = TRUE)
  )
  if (laplace) {
    out$insampleLogLikSubject <- insampleSubject
    out$outsampleLogLikSubject <- llrowoosSubject
    out$scoring <- if (subjectlevel) {
      "subject marginal at the fold's parameters; rows not scored"
    } else if (withheld) {
      "rows at random effect modes from the fold's training rows"
    } else "rows at random effect modes from the full data"
    if (corrected) {
      out$scoring <- paste0(out$scoring, "; Laplace optimum, not the fit's ",
        "quadrature-corrected estimate, and folds are not corrected")
    }
  }
  # Whichever route, a fold that was actually refit was refit with
  # `certify = FALSE` (see the call above): its estimate is not
  # curvature-certified, corrected, or checked against a better point the way
  # the full fit's is -- each fold would need its own Hessian, and folds exist
  # to be cheap. `scoring` is where a laplace fit already reports what a fold
  # does not get; an augmented fit's folds get the same caveat and this is the
  # only place a user reads it.
  foldstepped <- vapply(folded[usable], function(x) isTRUE(x$stepped), logical(1L))
  if (stepping) {
    out$foldinfluence <- vapply(folded[usable], function(x)
      as.numeric(x$influence)[1L], numeric(1L))
    out$foldstepped <- foldstepped
  }
  if (isTRUE(refit) && !all(foldstepped)) {
    uncertified <- paste0("folds are refit with optimcontrol$certify = FALSE ",
      "and are not curvature-certified or corrected")
    out$scoring <- if (is.null(out$scoring)) uncertified else
      paste0(out$scoring, "; ", uncertified)
  }
  if (any(foldstepped)) {
    stepped <- if (all(foldstepped)) {
      "folds are not refit: one Newton step from the full estimate"
    } else sprintf(paste0("%d of %d folds take one Newton step from the full ",
      "estimate (influence below %g), the rest are refit"), sum(foldstepped),
      length(foldstepped), casewiseBar)
    out$scoring <- if (is.null(out$scoring)) stepped else
      paste0(out$scoring, "; ", stepped)
  }
  out
}

# A refit fold stops once a step is predicted to gain less than this, in log
# likelihood: see the refit in `.ctBackendLOOFold()`.
.ctLOOFoldGain <- 1e-4

# `casewiseApproximation` as an influence bar: a fold whose step is shorter is
# stepped, the rest are refit. TRUE steps every fold, FALSE none.
.ctLOOCasewiseBar <- function(x) {
  if (isTRUE(x)) return(Inf)
  if (isFALSE(x)) return(0)
  if (is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 0) return(as.numeric(x))
  stop("casewiseApproximation must be TRUE, FALSE, or one non-negative ",
    "influence bar.", call. = FALSE)
}

# The score of a fold's held-out data at the estimate, without the prior: each
# unit the fold empties gives its whole likelihood score, and a unit it leaves
# partly observed the difference between its scores on the full and the
# training data. That costs a gradient of the training objective, spent only
# when the step may be the fold's answer (`training`, a function returning the
# fold's training fit, is given); otherwise NULL, as when there are no scores
# or the training scores cannot be had -- the fold is then refit from the
# estimate.
.ctLOOFoldScore <- function(ctx, heldout, training = NULL) {
  scores <- ctx$scores
  if (is.null(scores) || is.null(ctx$covariance)) return(NULL)
  rowunit <- ctx$unitOf[ctx$rowsubject]
  units <- unique(rowunit[heldout])
  emptied <- vapply(units, function(U) all(which(rowunit == U) %in% heldout),
    logical(1L))
  removed <- colSums(scores[units, , drop = FALSE])
  if (all(emptied)) return(removed)
  if (is.null(training)) return(NULL)
  trained <- try(.ctBackendScoreMatrix(training(), ctx$est, prior = FALSE),
    silent = TRUE)
  if (inherits(trained, "try-error")) return(NULL)
  removed - colSums(trained[units[!emptied], , drop = FALSE])
}

# One cross-validation fold: refit without the held-out rows (or take the
# score-based step instead), then score every row against the full data. Run in
# this session or in a worker process (`ctLOO(parallelFolds = TRUE)`), so it
# reads nothing but `job`: the held-out rows and the shared `context` built in
# `.ctBackendLOO()`. NULL when the refit failed.
#' @keywords internal
.ctBackendLOOFold <- function(job) {
  heldout <- job$heldout
  ctx <- job$context
  fit <- ctx$fit; spec <- ctx$spec; model <- ctx$model; est <- ctx$est
  rowsubject <- ctx$rowsubject; refit <- ctx$refit
  laplace <- ctx$laplace
  withheld <- ctx$withheld; subjectlevel <- ctx$subjectlevel
  covariance <- ctx$covariance; tol <- ctx$tol
  cores <- ctx$cores
  pars <- est
  training <- spec$data
  # `[-heldout, ]` would drop the rows; setting them missing keeps the
  # subject's timeline intact, so the filter still propagates across the
  # gap exactly as it will when the row is scored out of sample.
  training[heldout, model$manifestNames] <- NA
  # Built when something reads it: a refit, a step that needs the training
  # data's score, or a Laplace fit's training modes.
  trainingfit <- NULL
  trainingFit <- function() {
    if (is.null(trainingfit)) trainingfit <<- .ctFitReplaceData(fit, training)
    trainingfit
  }
  influence <- NA_real_
  stepped <- FALSE
  if (isTRUE(refit)) {
    start <- est
    removed <- .ctLOOFoldScore(ctx, heldout,
      training = if (isTRUE(ctx$stepping)) trainingFit else NULL)
    if (!is.null(removed)) {
      step <- as.numeric(covariance %*% removed)
      start <- est - step
      influence <- sum(removed * step)
    }
    stepped <- isTRUE(ctx$stepping) && is.finite(influence) &&
      influence < ctx$casewiseBar
    if (stepped) {
      pars <- start
    } else {
      # `tol` means the same thing here as it does on the stan branch above
      # and as `optimcontrol$tol` does in ctFit(): the objective tolerance.
      # It used to be handed to the engine's gradient criterion instead, so
      # the same argument relaxed two different things depending on backend.
      # Unless given it is not set. The fold stops instead once a step is
      # predicted to gain less than 1e-4 (`innergaptol`, in log likelihood):
      # a remaining gain G moves a held-out unit's term by at most about
      # sqrt(2 G h), h its influence, so 0.014 sqrt(h). Against converged
      # refits on a 40-subject random intercept with an outlier, that was
      # 1e-4 typically and 0.014 at worst, at the cost of the old objective
      # tolerance of 1e-5 (0.003 and 0.13), where converging took 23 times
      # as long; on a random-drift model 1e-5 was off by 1 in total.
      #
      # `certify = FALSE` is explicit rather than inherited, because a fold
      # is never certified or corrected after this call (see the comment on
      # the Laplace optimum in `.ctBackendLOO()`) -- that is deliberate, but
      # the fit's own
      # optimcontrol does not say so, so a fold used to inherit whatever
      # `certify` the full fit was made with. That reached
      # `.ctBackendInnerGapTol()`, which aims the optimiser's cheap stopping
      # rule inside the certification tolerance only because something will
      # spend a Hessian closing the rest of the gap; nothing does that for a
      # fold, so the rule is set for the fold's own purpose (above) rather
      # than inherited.
      # `saem` and `stochastic` off for the same kind of reason: both are
      # phases for a stage that starts away from the optimum, and a fold
      # starts at the full fit's estimate.
      result <- try(.ctJuliaOptimise(trainingFit()$model_spec, start,
        optimcontrol = utils::modifyList(
          as.list(fit$args$resolved$optimcontrol),
          c(if (!is.null(tol)) list(tol = tol),
            list(certify = FALSE, saem = FALSE, stochastic = FALSE,
              innergaptol = .ctLOOFoldGain))),
        cores = cores),
        silent = TRUE)
      if (inherits(result, "try-error")) return(NULL)
      pars <- as.numeric(result$minimizer)
    }
  }
  # Every row scored against the *full* data, so the held-out rows have a
  # likelihood to report and the rest are each fold's in-sample rows, as the
  # Stan path's are (it restores `dokalmanrows` before this step). When
  # `leaveOutN` is used without refitting, the in-sample likelihoods are the
  # ones wanted, and they are the same call -- one already made, when the
  # parameters are the estimate and no modes come from the training rows.
  #
  # Except after a step. Against a refit that pass is one evaluation among
  # many; against a step it is the fold's whole cost, and leaving each of N
  # subjects out made it N passes over N subjects -- an hour on 1046 subjects
  # where the step itself took a second. So a stepped fold on the filter's
  # route leaves its held-out rows to one pass shared by every stepped fold
  # (`.ctBackendLOOScoreSteps()`), and a stepped Laplace subject fold, whose
  # answer is its units' marginals below, scores no rows.
  deferred <- stepped && !laplace
  rowless <- stepped && subjectlevel
  modesFrom <- if (laplace && withheld && !rowless) trainingFit()$model_spec else NULL
  out <- list(llrow = if (deferred || rowless) rep(NA_real_, length(rowsubject)) else
    if (!isTRUE(refit) && is.null(modesFrom)) ctx$insample else
    .ctBackendRowLoglik(fit, pars, modesFrom = modesFrom), pars = pars,
    influence = influence, stepped = stepped, deferred = deferred)
  if (subjectlevel) {
    # The held-out units' marginals, and no other unit's: nothing reads them.
    held <- unique(rowsubject[heldout])
    full <- .ctBackendUnitTerms(spec, pars, subjects = held)
    train <- .ctBackendUnitTerms(trainingFit()$model_spec, pars,
      subjects = held)
    out$subject <- .ctBackendLOOHeldSubjects(full, train, held)
    out$converged <- isTRUE(all(full$converged)) && isTRUE(all(train$converged))
  }
  out
}

# The held-out rows of the folds a step answered on the filter's route
# (`deferred`, indices into `folded`), scored together: one filter pass with a
# parameter row per subject (`ctsem_kalman`'s second form), each subject at the
# parameters of the fold holding its rows. Subject folds hold each subject in
# one fold, so that is one pass for them all; a subject with rows in several
# folds (row folds) takes a pass for each.
.ctBackendLOOScoreSteps <- function(fit, est, folded, deferred, samplerows,
  rowsubject) {
  subjects <- unique(rowsubject)
  index <- match(rowsubject, subjects)
  pending <- lapply(deferred, function(f) unique(index[samplerows[[f]]]))
  while (any(lengths(pending))) {
    values <- matrix(est, length(subjects), length(est), byrow = TRUE)
    taken <- logical(length(subjects))
    scored <- vector("list", length(deferred))
    for (k in seq_along(deferred)) {
      s <- pending[[k]][!taken[pending[[k]]]]
      if (!length(s)) next
      values[s, ] <- matrix(folded[[deferred[k]]]$pars, length(s), length(est),
        byrow = TRUE)
      taken[s] <- TRUE
      scored[[k]] <- s
      pending[[k]] <- setdiff(pending[[k]], s)
    }
    llrow <- .ctBackendRowLoglik(fit, values)
    for (k in which(lengths(scored) > 0L)) {
      f <- deferred[k]
      rows <- samplerows[[f]][index[samplerows[[f]]] %in% scored[[k]]]
      folded[[f]]$llrow[rows] <- llrow[rows]
    }
  }
  folded
}


# Refit-free PSIS leave-one-out ----------------------------------------------
#
# Pareto-smoothed importance sampling, after bigIRT's `looIRT`. Two levels:
#
#   units  Leave one unit out: the outermost random-effect level, because
#          subjects sharing a group effect are not independent given theta,
#          and "leave one subject out" of one would need that subject's
#          marginal given its group mates, which the ratio trick cannot give.
#          A unit's term ll_i(theta) is its likelihood with its random effects
#          integrated out -- the filter's subject term on the augmented route,
#          the Laplace layer's unit term otherwise (`.ctBackendUnitTerms()`).
#          Population draws theta_s are a sampled fit's posterior draws, whose
#          ratios removing unit i are 1 / p(y_i | theta_s), or for an optimized
#          fit a proposal q = N(est, s^2 cov) corrected to the posterior by
#          lp_s - log q_s, so log r_is = lp_s - log q_s - ll_is.
#          elpd_i = log sum_s w_is exp(ll_is), w the smoothed, normalised
#          ratios. A unit whose k is above 0.7 moves the posterior further than
#          those draws reach, so it is estimated again from draws that include
#          N(est - cov s_i, s^2 cov) -- est - cov s_i is the
#          casewiseApproximation step, s_i the unit's likelihood score -- with
#          ratios to the posterior without it (`.ctPsisRefine()`). That needs lp
#          to be the posterior the draws describe, so not on a fit sampled
#          jointly over its random effects ('none'), whose unit terms are
#          Laplace approximations to the sampled marginal; nor when the
#          proposal itself misses the posterior (its own k above 0.7), which
#          makes every unit's k high and is a warning instead.
#   rows   Leave one row out at the estimate. Each unit's effects are drawn
#          from N(uhat_U, s^2 M_U^-1) and corrected to p(u | y) by
#          log v_s = g_U(u_s) - log q(u_s), g_U being the unit's rows plus its
#          standard normal prior; row r is removed by log r_s = log v_s -
#          llrow_rs. The llrow are one-step-ahead terms, so this is pointwise
#          LOO of the filter's factorisation, not exact time-series LOO.
#
# Both proposals are widened by `scale` = 1.5, and that is derived rather than
# tuned. The ratio for point i targets the posterior without i, whose precision
# is the full posterior's P less i's contribution c. With a proposal of
# precision P / scale^2 the ratios have finite variance only while
# 2 (P - c) > P / scale^2: at scale 1 that fails for any point carrying half
# the posterior's information, which a subject's first observation of a random
# intercept easily does. At 1.5 the bound is c < 0.78 P. Measured against the
# exact oracles in test-julia-loo-laplace.R (16 subjects, two free parameters,
# local machine, six seeds): at scale 1 the row-level rms error at 4000 draws
# was 0.010-0.019 with Pareto k up to 0.81; at 1.5 it was 0.004-0.008 with k
# below 0, about half its value at 1000 draws. Scale 2 was no better than 1.5
# at the unit level.

.ctLogSumExp <- function(x) {
  top <- max(x)
  if (!is.finite(top)) return(top)
  top + log(sum(exp(x - top)))
}

# Smoothed, self-normalised log weights and Pareto k for each column of a
# draws-by-points matrix of log ratios. `r_eff`, one per column, is the draws'
# relative efficiency (MCMC draws), 1 for independent ones.
.ctPsisWeights <- function(logratios, r_eff = NULL) {
  logratios <- as.matrix(logratios)
  if (is.null(r_eff)) r_eff <- rep(1, ncol(logratios))
  smoothed <- suppressWarnings(loo::psis(logratios, r_eff = r_eff))
  list(logweights = matrix(stats::weights(smoothed, log = TRUE, normalize = TRUE),
    nrow = nrow(logratios)), k = as.numeric(loo::pareto_k_values(smoothed)))
}

.ctBackendLOOPsis <- function(fit, subjectwise, ndraws, cores, thin = TRUE,
  scale = 1.5) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("method = 'psis' needs the loo package: install.packages('loo').",
      call. = FALSE)
  }
  ndraws <- suppressWarnings(as.integer(ndraws)[1L])
  if (is.na(ndraws) || ndraws < 50L) {
    stop("nsamples must be a whole number of at least 50.", call. = FALSE)
  }
  if (isFALSE(fit$args$resolved$intoverstates)) {
    stop("method = 'psis' needs a fit that integrates its states ",
      "(intoverstates = TRUE): a unit's term is its likelihood with its ",
      "states and random effects integrated out.", call. = FALSE)
  }
  spec <- .ctBackendSpec(fit)
  sampled <- !is.null(fit$sample)
  if (!isTRUE(subjectwise)) {
    if (!.ctBackendIsLaplace(fit) || sampled) {
      stop("method = 'psis' leaves one row out only on an optimized fit with ",
        "intoverpop = 'laplace', where each unit's random effects are drawn ",
        "from the Gaussian the Laplace approximation fits. Leave one subject ",
        "out (subjectwise = TRUE), or leave rows out by refitting with ",
        "method = 'kfold', subjectwise = FALSE.", call. = FALSE)
    }
    # The Laplace posterior is what the draws are corrected to, so they are
    # centred at its mode; see `.ctBackendLOO()` for the same choice.
    est <- .ctLaplaceOptimum(fit)
    return(.ctBackendWithMaxChunks(cores,
      .ctBackendLOOPsisRows(fit, spec, est, ndraws, scale)))
  }
  .ctBackendWithMaxChunks(cores,
    .ctBackendLOOPsisUnits(fit, spec, ndraws, scale, sampled, thin))
}

# Summaries shared by both levels, from pointwise elpd and lpd.
.ctLOOPsisSummary <- function(pointwise, level, ndraws, ndropped, proposal_k, note) {
  elpd <- pointwise$elpd_loo
  ok <- is.finite(elpd)
  n <- sum(ok)
  out <- list(
    elpd_loo = sum(elpd[ok]),
    se_elpd_loo = if (n > 1L) sqrt(n * stats::var(elpd[ok])) else NA_real_,
    p_loo = sum(pointwise$p_loo[ok]),
    looic = -2 * sum(elpd[ok]),
    pointwise = pointwise,
    n_high_k = sum(pointwise$pareto_k > 0.7, na.rm = TRUE),
    ndraws = ndraws, ndropped = ndropped, proposal_k = proposal_k,
    level = level, note = note)
  class(out) <- "ctLOOpsis"
  out
}

# Normal draws around `centre` with covariance `scale^2 covariance`, and their
# log density: the proposal the optimized fit's draws, and a refined unit's,
# come from.
.ctPsisProposal <- function(centre, covariance, ndraws, scale) {
  covariance <- as.matrix(covariance)
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  if (any(!is.finite(eig$values)) || any(eig$values <= 0)) {
    stop("The fit's covariance is not positive definite, so it cannot serve as ",
      "the proposal for method = 'psis'.", call. = FALSE)
  }
  npar <- length(centre)
  z <- matrix(stats::rnorm(npar * ndraws), npar, ndraws)
  list(draws = centre + scale * (eig$vectors %*% (sqrt(eig$values) * z)),
    logq = -npar / 2 * log(2 * pi) - sum(log(eig$values)) / 2 -
      npar * log(scale) - colSums(z^2) / 2)
}

# The chain of each of a sampled fit's posterior draws, in
# `fit$estimate$rawposterior`'s row order (chain after chain), or NULL where
# the counts do not add up to it.
.ctPsisChains <- function(fit, ndraws) {
  s <- fit$sample
  lengths <- if (!is.null(s$chain_draws)) as.integer(s$chain_draws) else
    rep(as.integer(s$draws)[1L], as.integer(s$chains)[1L])
  if (!length(lengths) || anyNA(lengths) || sum(lengths) != ndraws) return(NULL)
  rep(seq_along(lengths), lengths)
}

# Each column's relative efficiency over the chains (loo::relative_eff on the
# likelihood, which is scale free, so each column is shifted to its maximum
# first); 1 where it cannot be had, as for chains of unequal length.
.ctPsisRelEff <- function(ll, chain) {
  if (is.null(chain) || length(unique(chain)) < 2L ||
      length(unique(tabulate(chain))) > 1L) return(rep(1, ncol(ll)))
  lik <- exp(sweep(ll, 2L, apply(ll, 2L, max)))
  r <- try(suppressWarnings(loo::relative_eff(lik, chain_id = chain)), silent = TRUE)
  if (inherits(r, "try-error") || length(r) != ncol(ll) || any(!is.finite(r))) {
    return(rep(1, ncol(ll)))
  }
  as.numeric(r)
}

# Labels for the units, and the level's name: the outermost level's group
# labels under the Laplace layer, subject ids on the filter's route.
.ctPsisUnitLabels <- function(fit, spec, unit) {
  nunits <- max(unit)
  first <- match(seq_len(nunits), unit)
  idname <- .ctFitModelObject(fit)$subjectIDname
  ids <- unique(spec$data[[idname]])
  if (is.null(spec$laplace)) return(list(label = ids[first], name = idname))
  levels <- .ctSpecRandomEffectLevels(spec)
  outer <- levels[[length(levels)]]
  label <- if (!is.null(outer$labels) && length(outer$units) >= max(first)) {
    outer$labels[outer$units[first]]
  } else ids[first]
  name <- if (is.null(outer$name) || is.na(outer$name)) idname else outer$name
  list(label = label, name = name)
}

.ctBackendLOOPsisUnits <- function(fit, spec, ndraws, scale, sampled, thin) {
  if (sampled) {
    posterior <- as.matrix(fit$estimate$rawposterior)
    chain <- .ctPsisChains(fit, nrow(posterior))
    rows <- seq_len(nrow(posterior))
    # Evenly spaced through every chain, so each keeps its share.
    if (thin && ndraws < length(rows)) {
      rows <- unique(round(seq(1, length(rows), length.out = ndraws)))
    }
    draws <- t(posterior[rows, , drop = FALSE])
    if (!is.null(chain)) chain <- chain[rows]
    centre <- colMeans(posterior)
    covariance <- stats::cov(posterior)
    logq <- NULL
  } else {
    covariance <- fit$estimate$cov
    if (is.null(covariance)) {
      stop("method = 'psis' with subjectwise = TRUE draws from the fit's normal ",
        "approximation, and this fit carries no covariance (fit$estimate$cov). ",
        "Fit without optimcontrol = list(estonly = TRUE), or run ",
        "ctFitUncertainty(fit) first.", call. = FALSE)
    }
    covariance <- as.matrix(covariance)
    # The Laplace posterior is what the draws are corrected to, so they are
    # centred at its mode; see `.ctBackendLOO()` for the same choice.
    centre <- .ctLaplaceOptimum(fit)
    proposal <- .ctPsisProposal(centre, covariance, ndraws, scale)
    draws <- proposal$draws
    logq <- proposal$logq
    chain <- NULL
  }

  terms <- .ctBackendUnitTerms(spec, draws)
  keep <- is.finite(terms$value) & terms$converged &
    apply(terms$unit_loglik, 2L, function(x) all(is.finite(x)))
  ndropped <- sum(!keep)
  if (sum(keep) < 50L) {
    stop("Only ", sum(keep), " of ", ncol(draws), " draws could be evaluated; ",
      "method = 'psis' needs each unit's term at the population draws.",
      call. = FALSE)
  }
  if (ndropped) {
    warning(ndropped, " of ", ncol(draws), " draws were dropped because a ",
      "unit's term could not be evaluated there.", call. = FALSE)
  }

  ll <- t(terms$unit_loglik[, keep, drop = FALSE])      # draws x units
  nunits <- ncol(ll)
  if (sampled) {
    loo <- .ctPsisWeights(-ll, .ctPsisRelEff(ll, chain[keep]))
    lpdweights <- rep(-log(nrow(ll)), nrow(ll))
    proposal_k <- NA_real_
  } else {
    correction <- terms$value[keep] - logq[keep]         # q -> posterior
    loo <- .ctPsisWeights(correction - ll)
    posterior <- .ctPsisWeights(matrix(correction, ncol = 1L))
    lpdweights <- posterior$logweights[, 1L]
    proposal_k <- posterior$k
  }
  elpd <- vapply(seq_len(nunits), function(j)
    .ctLogSumExp(loo$logweights[, j] + ll[, j]), numeric(1L))
  lpd <- vapply(seq_len(nunits), function(j)
    .ctLogSumExp(lpdweights + ll[, j]), numeric(1L))
  k <- loo$k

  # A proposal that does not reach the posterior's tails makes every unit's k
  # high, whatever the unit does; redrawing units then fixes nothing, and the
  # fit's posterior is better sampled than approximated.
  mismatched <- isTRUE(proposal_k > 0.7)
  if (mismatched) {
    warning("The fit's normal approximation does not reach its posterior's ",
      "tails (Pareto k ", round(proposal_k, 2), "), so method = 'psis' is ",
      "unreliable for every unit. Sample the posterior -- ",
      "ctFitUncertainty(fit, 'sample') -- and run ctLOO(method = 'psis') on ",
      "that.", call. = FALSE)
  }
  refined <- rep(FALSE, nunits)
  refinable <- !sampled || identical(fit$sample$target, "marginal")
  high <- which(k > 0.7)
  if (length(high) && refinable && !mismatched) {
    base <- if (sampled) NULL else list(draws = draws[, keep, drop = FALSE],
      value = terms$value[keep], unit_loglik = terms$unit_loglik[, keep, drop = FALSE])
    again <- .ctPsisRefine(fit, spec, high, centre, covariance, ndraws, scale, base)
    better <- is.finite(again$elpd) & is.finite(again$k) & again$k < k[high]
    elpd[high[better]] <- again$elpd[better]
    k[high[better]] <- again$k[better]
    refined[high[better]] <- TRUE
  }

  labels <- .ctPsisUnitLabels(fit, spec, terms$unit)
  pointwise <- data.frame(unit = labels$label,
    nsubjects = tabulate(terms$unit, nbins = nunits),
    elpd_loo = elpd, lpd = lpd, p_loo = lpd - elpd, pareto_k = k,
    refined = refined, stringsAsFactors = FALSE)
  integrated <- if (is.null(spec$laplace)) "by the filter" else
    "by the Laplace approximation"
  note <- paste0("Leave one '", labels$name, "' out, its random effects ",
    "integrated ", integrated, "; population parameters ",
    if (sampled) "from the posterior draws." else paste0("integrated over the ",
      "fit's normal approximation, importance corrected to the posterior."))
  if (any(k > 0.7) && !refinable) {
    note <- paste0(note, " Units with Pareto k above 0.7 are not redrawn: the ",
      "fit sampled its random effects jointly, so its objective is not the ",
      "posterior its draws describe.")
  }
  .ctLOOPsisSummary(pointwise, level = labels$name, ndraws = sum(keep),
    ndropped = ndropped, proposal_k = proposal_k, note = note)
}

# Units whose Pareto k is too high, each estimated again by importance
# sampling to the posterior without it, lp - ll_i, from an equal mixture of
# two proposals with the same widened covariance: the one around `centre`
# (an optimized fit's own draws, `base`, reused; drawn here for a sampled fit),
# which covers the full posterior, and one shifted to where leaving the unit
# out moves the estimate -- one Newton step by the unit's likelihood score, as
# casewiseApproximation takes -- which covers the posterior without it. The
# ratio estimate is p(y | ...) / p(y_-i | ...), so it needs both. Ratios are to
# the mixture density (the balance heuristic). Measured on two-parameter
# fits with an outlying subject against a grid (8 seeds, 2000 draws a
# proposal), the outlier's rms error was 0.30 and 0.10 by plain PSIS, 0.16 and
# 0.08 from the shifted proposal alone, 0.12 and 0.05 from the mixture, which
# was no worse than plain PSIS for the other units. NA where the scores or the
# draws fail.
.ctPsisRefine <- function(fit, spec, units, centre, covariance, ndraws, scale,
  base = NULL) {
  out <- list(elpd = rep(NA_real_, length(units)), k = rep(NA_real_, length(units)))
  scores <- try(.ctBackendScoreMatrix(fit, centre, prior = FALSE), silent = TRUE)
  if (inherits(scores, "try-error")) return(out)
  if (is.null(base)) {
    drawn <- .ctPsisProposal(centre, covariance, ndraws, scale)
    terms <- .ctBackendUnitTerms(spec, drawn$draws)
    base <- list(draws = drawn$draws, value = terms$value,
      unit_loglik = terms$unit_loglik)
  }
  for (j in seq_along(units)) {
    U <- units[j]
    shifted <- centre - as.numeric(covariance %*% scores[U, ])
    proposal <- .ctPsisProposal(shifted, covariance, ndraws, scale)
    terms <- .ctBackendUnitTerms(spec, proposal$draws)
    draws <- cbind(base$draws, proposal$draws)
    value <- c(base$value, terms$value)
    llU <- c(base$unit_loglik[U, ], terms$unit_loglik[U, ])
    q0 <- .ctPsisLogDensity(draws, centre, covariance, scale)
    q1 <- .ctPsisLogDensity(draws, shifted, covariance, scale)
    logq <- log(0.5) + pmax(q0, q1) + log1p(exp(-abs(q0 - q1)))
    ok <- is.finite(value) & is.finite(llU)
    if (sum(ok) < 50L) next
    weights <- .ctPsisWeights(matrix(value[ok] - llU[ok] - logq[ok], ncol = 1L))
    out$elpd[j] <- .ctLogSumExp(weights$logweights[, 1L] + llU[ok])
    out$k[j] <- weights$k
  }
  out
}

# The log density of each column of `draws` under N(centre, scale^2 covariance).
.ctPsisLogDensity <- function(draws, centre, covariance, scale) {
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  z <- crossprod(eig$vectors, draws - centre) / (scale * sqrt(eig$values))
  -nrow(draws) / 2 * log(2 * pi) - sum(log(eig$values)) / 2 -
    nrow(draws) * log(scale) - colSums(z^2) / 2
}

.ctBackendLOOPsisRows <- function(fit, spec, est, ndraws, scale) {
  module <- .ctJuliaModule(spec$project)
  model <- .ctFitModelObject(fit)
  seed <- sample.int(.Machine$integer.max, 1L)
  res <- .ctBackendJuliaValue(module$ctsem_laplace_effect_draws(
    .ctJuliaObjective(fit, hybrid = TRUE), .ctJuliaNumericVector(est), as.integer(ndraws),
    seed = as.integer(seed), scale = as.numeric(scale)))
  llrow <- matrix(as.numeric(res$llrow), ncol = ndraws)      # rows x draws
  logq <- matrix(as.numeric(res$logq), ncol = ndraws)        # units x draws
  logprior <- matrix(as.numeric(res$logprior), ncol = ndraws)
  converged <- as.logical(res$converged)
  unitOfSubject <- as.integer(res$unit)
  rowsubject <- .ctFitRowSubject(fit)
  if (length(rowsubject) != nrow(llrow)) {
    stop("The filter reported ", nrow(llrow), " rows for ", length(rowsubject),
      " data rows.", call. = FALSE)
  }
  rowunit <- unitOfSubject[rowsubject]
  nunits <- nrow(logq)

  # g_U(u_s) - log q(u_s): the unit's rows and prior, against the proposal.
  unitll <- matrix(0, nunits, ndraws)
  summed <- rowsum(llrow, rowunit)
  unitll[as.integer(rownames(summed)), ] <- summed
  logv <- unitll + logprior - logq

  observed <- rowSums(!is.na(as.matrix(spec$data[, model$manifestNames,
    drop = FALSE]))) > 0L
  usable <- observed & converged[rowunit]
  if (any(observed & !converged[rowunit])) {
    warning(sum(!converged), " of ", nunits, " units' inner mode solve failed; ",
      "their ", sum(observed & !converged[rowunit]), " rows are left NA.",
      call. = FALSE)
  }
  rows <- which(usable)
  if (!length(rows)) stop("No row could be scored.", call. = FALSE)
  ll <- t(llrow[rows, , drop = FALSE])                        # draws x rows
  loo <- .ctPsisWeights(t(logv[rowunit[rows], , drop = FALSE]) - ll)
  units <- sort(unique(rowunit[rows]))
  posterior <- .ctPsisWeights(t(logv[units, , drop = FALSE]))
  column <- match(rowunit[rows], units)
  elpd <- vapply(seq_along(rows), function(j)
    .ctLogSumExp(loo$logweights[, j] + ll[, j]), numeric(1L))
  lpd <- vapply(seq_along(rows), function(j)
    .ctLogSumExp(posterior$logweights[, column[j]] + ll[, j]), numeric(1L))

  ids <- unique(spec$data[[model$subjectIDname]])
  pointwise <- data.frame(row = rows, subject = ids[rowsubject[rows]],
    time = as.numeric(spec$data[[model$timeName]][rows]),
    elpd_loo = elpd, lpd = lpd, p_loo = lpd - elpd, pareto_k = loo$k,
    stringsAsFactors = FALSE)
  proposal_k <- rep(NA_real_, nunits)
  proposal_k[units] <- posterior$k
  .ctLOOPsisSummary(pointwise, level = "row", ndraws = ndraws,
    ndropped = sum(!converged), proposal_k = proposal_k,
    note = paste0("Leave one row out, conditional on the estimated population ",
      "parameters, integrating over each unit's random effects. Rows are the ",
      "filter's one-step-ahead terms, so this is pointwise leave-one-out of ",
      "that factorisation, not exact time-series leave-one-out."))
}

#' @export
print.ctLOOpsis <- function(x, ...) {
  cat("PSIS leave-one-out, level:", x$level, "\n")
  cat(sprintf("  elpd_loo %.2f (se %.2f), p_loo %.2f, looic %.2f\n",
    x$elpd_loo, x$se_elpd_loo, x$p_loo, x$looic))
  cat(sprintf("  %d %s, %d draws", nrow(x$pointwise),
    if (identical(x$level, "row")) "rows" else "units", x$ndraws))
  if (x$ndropped) cat(",", x$ndropped, if (identical(x$level, "row"))
    "units dropped" else "draws dropped")
  cat("\n")
  cat("  Pareto k > 0.7:", x$n_high_k)
  refined <- sum(x$pointwise$refined, na.rm = TRUE)
  if (refined) cat(";", refined, "redrawn from the leave-one-out step")
  cat("\n")
  # How well the normal approximation the draws came from covers the
  # posterior: NA for a sampled fit, one per unit at the row level.
  proposal <- suppressWarnings(max(x$proposal_k, na.rm = TRUE))
  if (is.finite(proposal)) {
    cat(sprintf("  Pareto k of the normal approximation%s: %.2f%s\n",
      if (length(x$proposal_k) > 1L) " (largest unit)" else "", proposal,
      if (proposal > 0.7) ", too high: it misses the posterior's tails" else ""))
  }
  cat(" ", x$note, "\n")
  invisible(x)
}
