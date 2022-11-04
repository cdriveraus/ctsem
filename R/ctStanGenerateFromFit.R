#' Add a \code{$generated} object to ctstanfit object, with random data generated from posterior of ctstanfit object
#'
#' @param fit ctstanfit object
#' @param nsamples Positive integer specifying number of datasets to generate. 
#' @param fullposterior Logical indicating whether to sample from the full posterior (original nsamples) or the posterior mean.
#' @param verboseErrors if TRUE, print verbose output when errors in generation encountered.
#' @param cores Number of cpu cores to use.
#' @return Matrix of generated data -- one dataset per iteration, according to original time and missingness structure.
#' @export
#' @examples
#' if(w32chk()){
#'
#' gen <- ctStanGenerateFromFit(ctstantestfit, nsamples=3,fullposterior=TRUE,cores=1)
#' plot(gen$generated$Y[3,,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
#' }
ctStanGenerateFromFit<-function(fit,nsamples=200,fullposterior=FALSE, verboseErrors=FALSE,cores=2){
  
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object!')
  
  if(nsamples > ncol(fit$stanfit$rawposterior)) fit <- ctAddSamples(fit,nsamples = nsamples,cores=1)
  
  if(!fullposterior){
    umat=matrix(fit$stanfit$rawest,nrow=length(fit$stanfit$rawest),ncol=nsamples)
    } else umat=t(fit$stanfit$rawposterior)[,sample(1:nrow(fit$stanfit$rawposterior),size=nsamples),drop=FALSE]
  
  
  if(fit$setup$recompile) {
    message('Compilation needed -- compiling (usually ~ 1 min)')
    genm <- rstan::stan_model(model_code = 
        ctStanModelWriter(ctm = fit$ctstanmodel,
          gendata = TRUE,
          extratforms = fit$setup$extratforms,
          matsetup=fit$ctstanmodel$modelmats$matsetup))
  } else {
    genm <- stanmodels$ctsmgen
  }
  message('Generating data from ',ifelse(fullposterior,'posterior', 'posterior mean'))
  message('Using ',cores,'/', parallel::detectCores(),' logical CPU cores')
  standata <- fit$standata
  # standata$intoverstates=0L #why doesnt this work??
  standata$savescores <- 0L #have to disable for data generation in same structure as original
  # genf <- stan_reinitsf(genm,standata) 
  
  
  cs=suppressMessages(stan_constrainsamples(sm =genm,standata = standata,cores=cores,samples = t(umat),
    savescores = FALSE, savesubjectmatrices = FALSE,dokalman = TRUE, onlyfirstrow = FALSE,pcovn = FALSE))
  fit$generated$Y <- cs$Y #,c(2,1,3)) 
  fit$generated$llrow <- cs$llrow
  fit$generated$llrow[fit$generated$llrow==0]<-NA
  fit$generated$stanmodel <- genm
  
  dimnames( fit$generated$Y)<-list(
    sample=1:dim(fit$generated$Y)[1],
    row=1:dim(fit$generated$Y)[2],
    fit$ctstanmodel$manifestNames)
  
  fit$generated$Y[fit$generated$Y==99999] <- NA
  
  return(fit)
}
