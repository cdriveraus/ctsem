#' K fold cross validation for ctStanFit objects
#'
#' @param fit ctStanfit object
#' @param folds Number of cross validation splits to use -- 10 folds implies that the 
#' model is re-fit 10 times, each time to a data set with 1/10 of the observations randomly removed.
#' @param cores Number of processor cores to use. 
#' @param parallelFolds compute folds in parallel or use cores to finish single folds faster.
#' @param subjectwise drop random subjects instead of data rows?
#' @param keepfirstobs do not drop first observation (more stable estimates)
#' @return list
#' @export
#'
#' @examples
#' \donttest{ 
#' ctLOO(ctstantestfit())
#' }
ctLOO <- function(fit,folds=10,cores=2,parallelFolds=TRUE, 
  subjectwise=FALSE,keepfirstobs=TRUE){
  bootstrap <- FALSE
  if(!'ctStanFit' %in% class(fit)|| !'list' %in% class(fit$stanfit)) stop('Not an optimized ctStanFit object')
  
  message('Using ',cores,'/', parallel::detectCores(),' available CPU cores')
  
  samplerows <- which(c(1,diff(fit$standata$subject)) < ifelse(keepfirstobs,1,999)) #include all rows, or all except first obs per subject
  if(subjectwise) samplerows <- 1:fit$standata$nsubjects
  
  srows <- sample(samplerows,
    length(samplerows),# *ifelse(bootstrap,folds,1),
    replace=bootstrap)
  
  srows <- split(
    srows,
    sort(1:length(srows) %% folds))
  
  if(subjectwise) srows <- unlist(lapply(srows,function(x) which(fit$standata$subject %in% x)))
  
  sdat=fit$standata
  smodel <- fit$stanmodel
  init=fit$stanfit$rawest#rnorm(length(fit$stanfit$rawest),0,.01)#

  if(parallelFolds && cores > 1){
    clctsem <- parallel::makeCluster(spec = min(cores,folds))
    parallel::clusterExport(clctsem,c('sdat','smodel','init'),envir=environment())
    on.exit({parallel::stopCluster(clctsem)},add = TRUE)
  } else clctsem <- NA
  folded <- flexlapply(clctsem,X = srows,fn = function(x) {
    library(ctsem)
    
    sdat$dokalmanrowsdata[x] <- 0L
    sdat$dokalmanrowsata[-x] <- 1L
    # browser()
    e <- try(stanoptimis(standata = sdat,sm = smodel,init = init,
      # optimcontrol(list(stochastic=TRUE)),
      estonly = TRUE,cores=ifelse(parallelFolds,1,cores),plot=10,verbose=0))
    if('try-error' %in% class(e))   e <- try(stanoptimis(standata = sdat,sm = fit$stanmodel,init = fit$stanfit$rawest,
      estonly = TRUE,cores=ifelse(parallelFolds,1,cores),stochastic=FALSE))
    
    if('try-error' %in% class(e)) return(NA) else{
    sdat$savescores <- 1L
    sdat$dokalmanrowsdata <- fit$standata$dokalmanrowsdata
    smf <- stan_reinitsf(smodel,sdat)
    cp = rstan::constrain_pars(smf,e$rawest)
    return(list(llrow=cp$llrow,pars=e$rawest))
    }
  },cores = ifelse(parallelFolds,cores,1))
 # browser() 
  smf <- stan_reinitsf(smodel,sdat)
  cp = rstan::constrain_pars(smf,fit$stanfit$rawest)
  ee=unlist(lapply(1:folds,function(x) -sum(folded[[x]]$llrow)/fit$standata$ndatapoints))
  out <- list(
    foldrows=srows,
    foldpars = as.matrix(data.frame(lapply(folded,function(x) x$pars))),
    outsampleEntropyFolds = ee,
    insampleEntropy = -sum(cp$llrow)/fit$standata$ndatapoints,
    outsampleEntropyMean = mean(ee,na.rm=TRUE),
    outsampleEntropySD = sd(ee,na.rm=TRUE)
  )
  return( out)
}

