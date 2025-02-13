#' Fit and summarise a list of ctsem models 
#'
#' @param mlist Named list of models
#' @param datalong ctsem long format data
#' @param type 'ct' for continuous time or 'dt' for discrete time
#' @param cores number of cpu cores to use
#' @param summaryOutput Generate summary output into ctSummary folder? Large datasets can take some time. 
#' @param saveFits Save fit objects to working directory?
#' @param summaryArgs Additional arguments for ctSummarise.
#' @param prefix prefix for output files.
#' @param cv Perform k-fold cross validation?
#' @param cvArgs Additional arguments for ctLOO function used for cross validation.
#' @param ... Additional arguments for ctStanFit. 
#'
#' @return List containing a named list of model fits ($fits), and a compare object ($compare)
#' @export
#'
#' @examples
#' \dontrun{
#' sunspots<-data.frame(id=1,
#'   time=do.call(seq,(lapply(attributes(sunspot.year)$tsp,function(x) x))),
#'   sunspots=sunspot.year)
#' 
#'  ssmodel1 <- ctModel(type='omx', manifestNames='sunspots', Tpoints=3,
#'   latentNames=c('ss_level', 'ss_velocity'),
#'    LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
#'    DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
#'    MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
#'    MANIFESTVAR=diag(0,1), #As per original spec
#'    CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'    DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#'  
#'  ssmodel2 <- ssmodel1
#'  ssmodel2$LAMBDA[2] <- 0
#'  
#'  fits<-ctFitMultiModel(list(m1=ssmodel1,m2=ssmodel2),datalong = sunspots,
#'    summaryOutput = FALSE, saveFits = FALSE, cores=1, cv=TRUE, cvArgs=list(folds=5))
#'  print(fits$compare)
#' }
 
 
ctFitMultiModel <- function(mlist, datalong, prefix='',type='ct',cores=2, summaryOutput=TRUE, 
  saveFits = TRUE, summaryArgs = list(),cv=FALSE, cvArgs=list(),...){
  
  newfit <- function(model,name){ #function to convert old model to new form, fit with and without covariates, summarise, and save.
    if(class(model) %in% "ctsemInit") model <- ctStanModel(model,type = type) #convert to new model form  
    fit <- ctStanFit(datalong =datalong, ctstanmodel =model,cores=cores,...)
    
    if(summaryOutput){
      summaryArgs$cores <- cores
      summaryArgs$name <- paste0(prefix,name)
      summaryArgs$sf <- fit
      do.call(ctSummarise, summaryArgs)
    }
    if(saveFits) save(fit,file=paste0('fit_',name,'.rda'))
    return(fit)
  }
  
  
  mfit <- lapply(1:length(mlist),function(x) newfit(mlist[[x]],names(mlist)[x]))
  names(mfit) <- names(mlist)
  mcompare <- data.frame(t(sapply(mfit,function(x) c(
    npars=length(x$stanfit$rawest),
    loglik=x$stanfit$transformedparsfull$ll,
    aic=2*length(x$stanfit$rawest)-2*x$stanfit$transformedparsfull$ll,
    logprob=x$stanfit$optimfit$value))))
  
  # mcompare <- mcompare[order(mcompare$aic),] #disabled ordering - return in same order as models input
  
  if(summaryOutput){
  sink(file = paste0(prefix,'_compare.txt'))
  print(mcompare)
  sink()
  }
  
  if(cv){
    if(is.null(cvArgs$cores)) cvArgs$cores <- cores
    cv <- lapply(mfit,function(x){
      args <- cvArgs
      args$fit <- x
      do.call(ctLOO,args)
    })
    mcompare$OOSloglik <- lapply(cv,function(x) x$outsampleLogLik)
  }
  
  return(list(fits = mfit, compare=mcompare))
}
