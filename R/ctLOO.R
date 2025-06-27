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
#' @return list
#' @export
#'
#' @examples
#' \donttest{ 
#' ctLOO(ctstantestfit)
#' }
ctLOO <- function(fit, folds = 10, cores = 2, parallelFolds = FALSE, tol = 1e-5,
  subjectwise = ifelse(length(unique(fit$standata$subject)) >= folds, TRUE, FALSE),
  keepfirstobs = FALSE, leaveOutN = NA, refit = TRUE, casewiseApproximation = FALSE) {
  
  if (!'ctStanFit' %in% class(fit) || !length(fit$stanfit$stanfit@sim) == 0) 
    stop('Not an optimized ctStanFit object')
  if(is.na(as.integer(folds))) stop('Folds must be an integer')
  
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
    clctsem <- parallelly::makeClusterPSOCK(min(cores, folds))
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





