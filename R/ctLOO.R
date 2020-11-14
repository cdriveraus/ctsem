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
#' @return list
#' @export
#'
#' @examples
#' \donttest{ 
#' ctLOO(ctstantestfit)
#' }
ctLOO <- function(fit,folds=10,cores=2,parallelFolds=FALSE, 
  subjectwise=ifelse(length(unique(fit$standata$subject)) > folds, TRUE, FALSE),
  keepfirstobs=FALSE, init='fit'){
  bootstrap <- FALSE
  if(!'ctStanFit' %in% class(fit)|| !length(fit$stanfit$stanfit@sim)==0) stop('Not an optimized ctStanFit object')
  
  message('Using ',cores,'/', parallel::detectCores(),' available CPU cores')
  
  samplerows <- which(c(1,diff(fit$standata$subject)) < ifelse(keepfirstobs,1,999)) #include all rows, or all except first obs per subject
  if(subjectwise) samplerows <- 1:fit$standata$nsubjects
  
  srows <- sample(samplerows,
    length(samplerows),# *ifelse(bootstrap,folds,1),
    replace=bootstrap)
  
  srows <- split(
    srows,
    sort(1:length(srows) %% folds))
  
  if(subjectwise) srows <- lapply(srows,function(x) which(fit$standata$subject %in% x))
  
  sdat=fit$standata
  smodel <- fit$stanmodel
  if(init=='fit') init=fit$stanfit$rawest else init=rnorm(length(fit$stanfit$rawest),0,.01)#

  if(parallelFolds && cores > 1){
    clctsem <- parallel::makeCluster(spec = min(cores,folds))
    parallel::clusterExport(clctsem,c('sdat','smodel','init','parallelFolds'),envir=environment())
    on.exit({parallel::stopCluster(clctsem)},add = TRUE)
  } else clctsem <- NA

  folded <- flexlapply(clctsem,X = srows,fn = function(x) {
    require(ctsem)
    sdat$dokalmanrows[x] <- 0L
    sdat$dokalmanrows[-x] <- 1L
    e <- try(stanoptimis(standata = sdat,sm = smodel,init = init,
      # optimcontrol(list(stochastic=TRUE)),
      estonly = TRUE,cores=ifelse(parallelFolds,1,cores),plot=10,verbose=0))
    # 
    # if('try-error' %in% class(e) && init='fit')   e <- try(stanoptimis(standata = sdat,sm = smodel,
    #   init = 0,
    #   estonly = TRUE,cores=ifelse(parallelFolds,1,cores),stochastic=FALSE))
    
    if('try-error' %in% class(e)) return(NA) else{
    sdat$savescores <- 1L
    sdat$dokalmanrows <- fit$standata$dokalmanrows
    smf <- stan_reinitsf(smodel,sdat)
    cp = rstan::constrain_pars(smf,e$rawest)
    return(list(llrow=cp$llrow,pars=e$rawest))
    }
  },cores = ifelse(parallelFolds,cores,1))
 # browser()
 # sdat$savescores <- 1L
 #  smf <- stan_reinitsf(smodel,sdat)
  # cp = rstan::constrain_pars(smf,fit$stanfit$rawest)
  llrowoos=unlist(lapply(1:folds,function(x) folded[[x]]$llrow[srows[[x]] ]))
  llrowoos=llrowoos[match(1:(fit$standata$ndatapoints),unlist(srows))]
  llrowoosSubject=sapply(unique(fit$standata$subject),function(x) 
    sum(llrowoos[fit$standata$subject==x],na.rm=FALSE) )
  llrow=fit$stanfit$transformedparsfull$llrow
  llrowSubject=sapply(unique(fit$standata$subject),function(x) 
    sum(llrow[fit$standata$subject==x],na.rm=TRUE) )
  plot(llrow,llrowoos,col=fit$standata$subject,pch=16)
  abline(b = 1,a=0)
  plot(density(llrow-llrowoos),main='Original - OOS LogLik Difference')
  abline(v=mean(llrow-llrowoos))
  # ee=unlist(lapply(1:folds,function(x) -sum(folded[[x]]$llrow)/fit$standata$ndatapoints))
  out <- list(
    foldrows=srows,
    foldpars = as.matrix(data.frame(lapply(folded,function(x) x$pars))),
    # outsampleLogLikFolds=lloos,
    insampleLogLikRow=llrow,
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
  return( out)
}

