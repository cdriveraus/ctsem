#' Add a \code{$generated} object to ctstanfit object, with random data generated from posterior of ctstanfit object
#'
#' Add a \code{$generated} object to ctstanfit object, with random data
#' generated from posterior of ctstanfit object.
#' \code{ctStanGenerateFromFit} is maintained as a backward-compatible alias.
#'
#' @param fit ctstanfit object
#' @param nsamples Positive integer specifying number of datasets to generate. 
#' @param fullposterior Logical indicating whether to sample from the full posterior (original nsamples) or the posterior mean.
#' @param verboseErrors if TRUE, print verbose output when errors in generation encountered.
#' @param cores Number of cpu cores to use.
#' @param intoverstates For \code{backend='julia'} fits: which route generates
#' the latent states. \code{'fit'} (the default) mirrors the fit, so a
#' state-explicit fit gets a state-explicit predictive and every other fit
#' generates through the filter's one-step-ahead predictive. That keeps the
#' identity that a generated dataset's \code{llrow} is the likelihood reported
#' while generating it.
#'
#' \code{FALSE} resamples the trajectory whatever the fit did: drawn from the
#' process at each subject's own parameters, then each observation from its
#' conditional distribution given the state at its row. That is a draw from the
#' model rather than from the filter's approximation of it, and the two differ
#' for a non-Gaussian indicator -- where the measurement update is an
#' assumed-density projection -- or for state-dependent dynamics.
#' \code{TRUE} forces the filter route.
#'
#' Either way the states are regenerated rather than reused: what a fit
#' contributes is its sampled individual differences, and the trajectory is
#' drawn afresh conditional on the subject-specific model.
#' @return Matrix of generated data -- one dataset per iteration, according to original time and missingness structure.
#' @aliases ctStanGenerateFromFit
#' @export
#' @examples
#' gen <- ctGenerateFromFit(ctstantestfit, nsamples=3,fullposterior=TRUE,cores=1)
#' plot(gen$generated$Y[3,,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
ctGenerateFromFit<-function(fit,nsamples=200,fullposterior=FALSE, verboseErrors=FALSE,cores=2,
  intoverstates='fit'){
  
  # The julia engine generates from the same forward pass it filters with,
  # drawing each row from its own prior predictive (see R/ctBackendKalman.R and
  # the engine's kalman_trace.jl).
  if(inherits(fit,'ctJuliaFit')){
    return(.ctBackendGenerateFromFit(fit,nsamples=nsamples,
      fullposterior=fullposterior,cores=cores,intoverstates=intoverstates))
  }
  if(!'ctStanFit' %in% class(fit)) stop('Not a ctStanFit object!')
  # Refused by name rather than accepted and ignored: the generated Stan model
  # has one generation route, so there is no choice here to honour, and an
  # argument that means something on one backend and nothing on the other is
  # the shape of mistake this package keeps paying for.
  if(!identical(as.character(intoverstates)[1],'fit')) stop(
    "intoverstates is a backend='julia' option: the generated Stan model ",
    'generates through its own filter and offers no alternative route. Drop ',
    'the argument, or refit with backend="julia".', call.=FALSE)
  
  # `nrow`, not `ncol`: `rawposterior` is samples by parameters, so both of
  # these asked whether more datasets were wanted than the model has
  # parameters. And the top-up was gated on `is.null(fit$stanfit$stanfit)`,
  # which `stanoptimis()` never leaves NULL -- it always stores a reinitialised
  # stan model object there -- so the branch could not fire and a fit with
  # twenty draws generated two hundred datasets by resampling them with
  # replacement, silently, whenever the parameter count happened to be the
  # smaller number.
  #
  # Optimized fits only. A sampled fit's draws are the posterior; there is no
  # covariance to draw more of them from, and asking for more datasets than
  # draws resamples them with replacement below, as it always did.
  if(fullposterior && nsamples > nrow(fit$stanfit$rawposterior) &&
      length(fit$stanfit$stanfit@sim) == 0) {
    fit <- suppressMessages(ctOptimUncertainty(fit, uncertainty='stored',
      finishsamples=nsamples, cores=1))
  }

  #if nsamples still larger than available, use replacement
  replace <- nsamples > nrow(fit$stanfit$rawposterior)

  if(!fullposterior){
    umat=matrix(fit$stanfit$rawest,nrow=length(fit$stanfit$rawest),ncol=nsamples)
    } else umat=t(fit$stanfit$rawposterior)[,sample(1:nrow(fit$stanfit$rawposterior),size=nsamples,replace = replace),drop=FALSE]
  
  
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

#' @export
ctStanGenerateFromFit <- ctGenerateFromFit
