#' Cross validation for ctsem fits: K fold by refitting, or PSIS for Laplace fits
#'
#' @param fit ctStanfit object
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
#' @param casewiseApproximation if TRUE, use a bootstrapped gradient contributions approach to approximate the cross validation parameters -- much faster but less reliable.
#' @param tol tolerance for optimisation of refitted samples, can generally be more relaxed than the tolerance used for fitting initially.
#' @param method \code{'kfold'} (the default) refits the model once per fold.
#'   \code{'psis'} is refit-free Pareto-smoothed importance sampling
#'   leave-one-out, available only for \code{backend = 'julia'} fits with
#'   \code{intoverpop = 'laplace'}, and refused for any other fit; it needs the
#'   \pkg{loo} package. The K-fold arguments \code{folds}, \code{parallelFolds},
#'   \code{tol}, \code{keepfirstobs}, \code{leaveOutN}, \code{refit} and
#'   \code{casewiseApproximation} mean nothing to it and are refused if given.
#' @param ndraws Number of importance sampling draws for \code{method = 'psis'}.
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
#' \strong{\code{method = 'psis'}.} Pareto-smoothed importance sampling
#' (Vehtari, Gelman and Gabry, 2017) over draws built from the Laplace
#' approximation, without refitting, at one of two levels:
#' \itemize{
#'   \item \code{subjectwise = TRUE}: leave one unit out, where a unit is the
#'   outermost level of random effect -- a subject in a one-level model, a
#'   group in a multilevel one, since subjects sharing a group effect are not
#'   independent given the population parameters. Population parameter draws
#'   come from the fit's normal approximation (\code{fit$estimate$cov}, so the
#'   fit must carry one) widened by a factor of 1.5, are importance corrected
#'   to the Laplace posterior, and are weighted by each unit's Laplace marginal
#'   to remove it.
#'   \item \code{subjectwise = FALSE}: leave one row out, conditional on the
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
#' @return For \code{method = 'kfold'}, a list with \code{foldrows},
#'   \code{foldpars}, \code{insampleLogLikRow}, \code{LogLikRowFolds},
#'   \code{outsampleLogLikRow}, \code{insampleLogLik}, \code{outsampleLogLik}
#'   and entropy and standard deviation summaries. An \code{intoverpop =
#'   'laplace'} fit adds \code{insampleLogLikSubject},
#'   \code{outsampleLogLikSubject} and \code{scoring}, naming what was scored
#'   (see Details). For \code{method = 'psis'}, an object of class
#'   \code{ctLOOpsis}: \code{elpd_loo}, \code{se_elpd_loo}, \code{p_loo},
#'   \code{looic}, \code{pointwise} (a data frame of \code{elpd_loo},
#'   \code{lpd}, \code{p_loo} and \code{pareto_k} per unit or row),
#'   \code{n_high_k}, \code{ndraws}, \code{ndropped}, \code{proposal_k} (the
#'   Pareto k of the correction from proposal to Laplace posterior: one value at
#'   the unit level, one per unit at the row level), \code{level} and
#'   \code{note}.
#' @references Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian
#'   model evaluation using leave-one-out cross-validation and WAIC.
#'   \emph{Statistics and Computing}, 27, 1413--1432.
#' @export
#'
#' @examples
#' \donttest{
#' ctLOO(ctstantestfit)
#' }
ctLOO <- function(fit, folds = 10, cores = 2, parallelFolds = FALSE, tol = 1e-5,
  subjectwise = ifelse(length(unique(.ctFitRowSubject(fit))) >= folds, TRUE, FALSE),
  keepfirstobs = FALSE, leaveOutN = NA, refit = TRUE, casewiseApproximation = FALSE,
  method = c("kfold", "psis"), ndraws = 500) {

  method <- match.arg(method)
  if (identical(method, "psis")) {
    if (!.ctBackendIsLaplace(fit)) {
      stop("method = 'psis' needs a backend = 'julia' fit with intoverpop = ",
        "'laplace'; use method = 'kfold' for this fit.", call. = FALSE)
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
    return(.ctBackendLOOPsis(fit, subjectwise = subjectwise, ndraws = ndraws,
      cores = cores))
  }
  if (!missing(ndraws)) {
    stop("ndraws is used only by method = 'psis'.", call. = FALSE)
  }

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

# Whether a fit integrates its random effects by Laplace. A fit that samples
# them (`intoverpop = 'none'`) carries the same `spec$laplace` description, so
# the structure's presence is not enough; see `.ctBackendIntOverPop`.
.ctBackendIsLaplace <- function(fit) {
  if (!inherits(fit, "ctJuliaFit") || !.ctFitIsJulia(fit)) return(FALSE)
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

# The Laplace objective at each column of `draws` (`npar x ndraws`), with each
# unit's own term. A unit is the outermost level of random effect: a subject in
# a one-level model, a group in a nested one. One bridge call for all columns.
.ctBackendLaplaceUnitTerms <- function(spec, draws) {
  draws <- as.matrix(draws)
  storage.mode(draws) <- "double"
  module <- .ctJuliaModule(spec$project)
  res <- .ctBackendJuliaValue(module$ctsem_laplace_unit_terms(
    .ctJuliaObjective(.ctBackendAsModel(spec)), JuliaConnectoR::juliaPut(draws)))
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

  # A Laplace fit's random effects are integrated out, so the fold's training
  # rows have to be kept away from them as well as from the parameters. The
  # condition matches the stan path's: the rows are withheld unless `leaveOutN`
  # is used without refitting, where the in-sample likelihoods are the ones
  # wanted. Subject folds on such a fit score each held-out subject's marginal
  # rather than its rows; see the Details of ?ctLOO.
  laplace <- .ctBackendIsLaplace(fit)
  withheld <- all(is.na(leaveOutN)) || isTRUE(refit)
  subjectlevel <- laplace && isTRUE(subjectwise) && all(is.na(leaveOutN))

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
    training <- spec$data
    # `[-heldout, ]` would drop the rows; setting them missing keeps the
    # subject's timeline intact, so the filter still propagates across the
    # gap exactly as it will when the row is scored out of sample.
    training[heldout, model$manifestNames] <- NA
    refitting <- isTRUE(refit) && !isTRUE(casewiseApproximation)
    trainingfit <- if (refitting || (laplace && withheld)) {
      .ctFitReplaceData(fit, training)
    } else NULL
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
        # `tol` means the same thing here as it does on the stan branch above
        # and as `optimcontrol$tol` does in ctFit(): the objective tolerance.
        # It used to be handed to the engine's gradient criterion instead, so
        # the same argument relaxed two different things depending on backend.
        result <- try(.ctJuliaOptimise(trainingfit$model_spec, start,
          optimcontrol = utils::modifyList(
            as.list(fit$args$resolved$optimcontrol), list(tol = tol)),
          cores = cores),
          silent = TRUE)
        if (inherits(result, "try-error")) return(NULL)
        pars <- as.numeric(result$minimizer)
      }
    }
    # Scored against the *full* data, so the held-out rows have a likelihood to
    # report; the Stan path restores `dokalmanrows` before this step for the
    # same reason. When `leaveOutN` is used without refitting, the in-sample
    # likelihoods are the ones wanted, and they are the same call.
    modesFrom <- if (laplace && withheld) trainingfit$model_spec else NULL
    out <- list(llrow = .ctBackendRowLoglik(fit, pars, modesFrom = modesFrom),
      pars = pars)
    if (subjectlevel) {
      full <- .ctBackendLaplaceUnitTerms(spec, pars)
      train <- .ctBackendLaplaceUnitTerms(trainingfit$model_spec, pars)
      out$subject <- .ctBackendLOOHeldSubjects(full, train,
        unique(rowsubject[heldout]))
      out$converged <- isTRUE(all(full$converged)) && isTRUE(all(train$converged))
    }
    out
  })

  usable <- !vapply(folded, is.null, logical(1L))
  if (!any(usable)) stop("Every cross-validation fold failed to refit.", call. = FALSE)
  if (any(!usable)) {
    warning(sum(!usable), " of ", folds, " folds failed to refit and are dropped.",
      call. = FALSE)
  }

  llrow <- .ctBackendRowLoglik(fit, est)
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
    fitted <- .ctBackendLaplaceUnitTerms(spec, est)
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
  }
  out
}


# Refit-free PSIS leave-one-out for Laplace fits -----------------------------
#
# Pareto-smoothed importance sampling over draws the Laplace approximation
# already describes, after bigIRT's `looIRT`. Two levels:
#
#   units  Leave one unit out. Population draws theta_s ~ q = N(est, s^2 cov) are
#          corrected to the Laplace posterior, lp_s - log q_s, and unit i is
#          removed by dividing by its marginal: log r_is = lp_s - log q_s -
#          ll_is. elpd_i = log sum_s w_is exp(ll_is) with w the smoothed,
#          normalised ratios. A unit is the outermost random-effect level,
#          because subjects sharing a group effect are not independent given
#          theta, and "leave one subject out" of one would need that subject's
#          marginal given its group mates, which the ratio trick cannot give.
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
# draws-by-points matrix of log ratios.
.ctPsisWeights <- function(logratios) {
  logratios <- as.matrix(logratios)
  smoothed <- suppressWarnings(loo::psis(logratios,
    r_eff = rep(1, ncol(logratios))))
  list(logweights = matrix(stats::weights(smoothed, log = TRUE, normalize = TRUE),
    nrow = nrow(logratios)), k = as.numeric(loo::pareto_k_values(smoothed)))
}

.ctBackendLOOPsis <- function(fit, subjectwise, ndraws, cores, scale = 1.5) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("method = 'psis' needs the loo package: install.packages('loo').",
      call. = FALSE)
  }
  ndraws <- suppressWarnings(as.integer(ndraws)[1L])
  if (is.na(ndraws) || ndraws < 50L) {
    stop("ndraws must be a whole number of at least 50.", call. = FALSE)
  }
  spec <- .ctBackendSpec(fit)
  est <- as.numeric(fit$estimate$raw)
  .ctBackendWithMaxChunks(cores, if (isTRUE(subjectwise)) {
    .ctBackendLOOPsisUnits(fit, spec, est, ndraws, scale)
  } else {
    .ctBackendLOOPsisRows(fit, spec, est, ndraws, scale)
  })
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

.ctBackendLOOPsisUnits <- function(fit, spec, est, ndraws, scale) {
  covariance <- fit$estimate$cov
  if (is.null(covariance)) {
    stop("method = 'psis' with subjectwise = TRUE draws from the fit's normal ",
      "approximation, and this fit carries no covariance (fit$estimate$cov). ",
      "Fit without optimcontrol = list(estonly = TRUE), or run ",
      "ctOptimUncertainty(fit) first.", call. = FALSE)
  }
  covariance <- as.matrix(covariance)
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  if (any(!is.finite(eig$values)) || any(eig$values <= 0)) {
    stop("The fit's covariance is not positive definite, so it cannot serve as ",
      "the proposal for method = 'psis'.", call. = FALSE)
  }
  npar <- length(est)
  z <- matrix(stats::rnorm(npar * ndraws), npar, ndraws)
  draws <- est + scale * (eig$vectors %*% (sqrt(eig$values) * z))
  logq <- -npar / 2 * log(2 * pi) - sum(log(eig$values)) / 2 - npar * log(scale) -
    colSums(z^2) / 2

  terms <- .ctBackendLaplaceUnitTerms(spec, draws)
  keep <- is.finite(terms$value) & terms$converged &
    apply(terms$unit_loglik, 2L, function(x) all(is.finite(x)))
  ndropped <- sum(!keep)
  if (sum(keep) < 50L) {
    stop("Only ", sum(keep), " of ", ndraws, " draws could be evaluated; ",
      "method = 'psis' needs the Laplace inner solve to succeed across the ",
      "fit's normal approximation.", call. = FALSE)
  }
  if (ndropped) {
    warning(ndropped, " of ", ndraws, " draws were dropped because the inner ",
      "mode solve failed there.", call. = FALSE)
  }

  ll <- t(terms$unit_loglik[, keep, drop = FALSE])      # draws x units
  correction <- terms$value[keep] - logq[keep]           # q -> Laplace posterior
  loo <- .ctPsisWeights(correction - ll)
  posterior <- .ctPsisWeights(matrix(correction, ncol = 1L))
  nunits <- ncol(ll)
  elpd <- vapply(seq_len(nunits), function(j)
    .ctLogSumExp(loo$logweights[, j] + ll[, j]), numeric(1L))
  lpd <- vapply(seq_len(nunits), function(j)
    .ctLogSumExp(posterior$logweights[, 1L] + ll[, j]), numeric(1L))

  levels <- .ctSpecRandomEffectLevels(spec)
  outer <- levels[[length(levels)]]
  first <- match(seq_len(nunits), terms$unit)
  idname <- .ctFitModelObject(fit)$subjectIDname
  label <- if (!is.null(outer$labels) && length(outer$units) >= max(first)) {
    outer$labels[outer$units[first]]
  } else unique(spec$data[[idname]])[first]
  pointwise <- data.frame(unit = label,
    nsubjects = tabulate(terms$unit, nbins = nunits),
    elpd_loo = elpd, lpd = lpd, p_loo = lpd - elpd, pareto_k = loo$k,
    stringsAsFactors = FALSE)
  unitname <- if (is.null(outer$name) || is.na(outer$name)) idname else outer$name
  .ctLOOPsisSummary(pointwise, level = unitname, ndraws = sum(keep),
    ndropped = ndropped, proposal_k = posterior$k,
    note = paste0("Leave one '", unitname, "' out: population parameters ",
      "integrated over the fit's normal approximation, importance corrected ",
      "to the Laplace posterior."))
}

.ctBackendLOOPsisRows <- function(fit, spec, est, ndraws, scale) {
  module <- .ctJuliaModule(spec$project)
  model <- .ctFitModelObject(fit)
  seed <- sample.int(.Machine$integer.max, 1L)
  res <- .ctBackendJuliaValue(module$ctsem_laplace_effect_draws(
    .ctJuliaObjective(fit), .ctJuliaNumericVector(est), as.integer(ndraws),
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
  cat("  Pareto k > 0.7:", x$n_high_k, "\n")
  cat(" ", x$note, "\n")
  invisible(x)
}
