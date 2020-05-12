#' Add a \code{$generated} object to ctstanfit object, with random data generated from posterior of ctstanfit object
#'
#' @param fit ctstanfit object
#' @param nsamples Positive integer specifying number of datasets to generate. 
#' @param fullposterior Logical indicating whether to sample from the full posterior (original nsamples) or the posterior mean.
#' @param verboseErrors if TRUE, print verbose output when errors in generation encountered.
#' @return Matrix of generated data -- one dataset per iteration, according to original time and missingness structure.
#' @export
#' @examples
#' \donttest{
#' if (!exists("ctstantestfit")) ctstantestfit <- ctstantestfitgen()
#' gen <- ctStanGenerateFromFit(ctstantestfit, nsamples=3,fullposterior=TRUE)
#' plot(gen$generated$Y[,3,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
#' }
ctStanGenerateFromFit<-function(fit,nsamples=200,fullposterior=FALSE, verboseErrors=FALSE){
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object!')
  if(!'stanfit' %in% class(fit$stanfit)) {
    umat=t(fit$stanfit$rawposterior)
  } else  {
    umat <- stan_unconstrainsamples(fit$stanfit,fit$standata)
  }

  
  if(!fullposterior) umat=matrix(apply(umat, 1, mean),ncol=1)
  umat=umat[,sample(1:ncol(umat),nsamples,replace = ifelse(nsamples > ncol(umat), TRUE, FALSE)),drop=FALSE]
  
  if(fit$setup$recompile) {
    message('Compilation needed -- compiling (usually ~ 1 min)')
    genm <- rstan::stan_model(model_code = 
        ctStanModelWriter(ctm = fit$ctstanmodel,
          gendata = TRUE,
          extratforms = fit$setup$extratforms,
          matsetup=fit$setup$matsetup))
  } else {
    genm <- stanmodels$ctsmgen
  }
  message('Generating data from posterior')
  standata <- fit$standata
  standata$savescores <- 0L #have to disable for data generation in same structure as original
  genf <- stan_reinitsf(genm,standata) 
  
  fit$generated$Y <- array(apply(umat, 2, function(x){
    out <- try(rstan::constrain_pars(genf,x)$Y)
    if('try-error' %in% class(out)) {
      out <- rep(NA, standata$ndatapoints)
      if(verboseErrors) {
        standata$verbose <<- 2L
        genf <- stan_reinitsf(genm,standata) 
        rstan::constrain_pars(genf,x)$Y
      }
    }
    return(out)
    }),dim=c(nrow(fit$standata$Y),fit$ctstanmodel$n.manifest,ncol(umat)))
  
  fit$generated$Y <- aperm(fit$generated$Y,c(1,3,2))
  fit$generated$Y[fit$generated$Y==99999] <- NA
  # print(fit$generated$Y)
  # browser()
  # for(i in 1:dim(fit$generated$Y)[3]){ #remove crazy outliers
  #   for(j in 1:dim(fit$generated$Y)[2]){
  #     # browser()
  #     fit$generated$Y[(fit$generated$Y[,j,i]) > 
  #         quantile(c(fit$generated$Y[,,i]),.99,na.rm=TRUE) * 100,j,i] <- NA
  #     fit$generated$Y[(fit$generated$Y[,j,i] < 
  #         quantile(c(fit$generated$Y[,,i]),.01,na.rm=TRUE) ),j,i] <- NA
  #   }
  # }
  return(fit)
}
