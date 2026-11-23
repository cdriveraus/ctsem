#' Generate data from a fitted model
#'
#' Adds a \code{$generated} object to a fit, holding random data generated
#' from the fit's posterior (or its estimate), on either backend.
#' \code{ctStanGenerateFromFit} is maintained as a backward-compatible alias.
#'
#' @param fit A fit from \code{\link{ctFit}}, on either backend.
#' @param nsamples Positive integer specifying number of datasets to generate. 
#' @param fullposterior Logical indicating whether to sample from the full posterior (original nsamples) or the posterior mean.
#' @param verboseErrors if TRUE, print verbose output when errors in generation encountered.
#' @param cores Number of cpu cores to use: R worker processes for a stan fit,
#' and for a julia fit the engine threads its draws are generated on.
#' @param intoverstates For \code{backend='julia'} fits, how the latent states
#' are generated. \code{'fit'} (the default) follows the fit. A fit with
#' explicit states (\code{intoverstates=FALSE}) then generates states
#' explicitly, and every other fit generates each observation from the Kalman
#' filter's one-step-ahead prediction, the density the fit used. A generated
#' dataset's \code{llrow} is then the likelihood under which it was
#' generated.
#'
#' \code{FALSE} simulates the model directly, whatever the fit used: each
#' subject's trajectory from the process at its own parameters, and then each
#' observation from its measurement model given the state, as
#' \code{\link{ctGenerate}} does. The two differ for non-Gaussian indicators
#' and for state dependent dynamics, where the filter is an approximation.
#' Comparing data generated each way, with the same \code{effects}, shows the
#' effect of that approximation, separately from any misfit of the model to
#' the observed data. \code{TRUE} always generates through the filter.
#'
#' In all cases the states are generated anew.
#' @param effects For \code{backend='julia'} fits with random effects.
#' \code{'population'} (the default) draws new subjects, and groups at every
#' level, from the fitted population distribution at each parameter draw, so
#' the generated data are the model's prediction for a new sample with the
#' same design. \code{'fitted'} keeps the fitted subjects, each at its sampled
#' effects when the fit sampled and saved them, and otherwise at its
#' conditional mode. A fit with \code{intoverpop='augmented'} has no
#' per-subject effects, so only \code{'population'} applies. Stan fits
#' generate as the fit was made, and give an error if the argument is set.
#' @return The fit, with a \code{$generated} element: \code{Y}, the generated
#' data as an array of sample by data row by manifest variable, on the
#' original time and missingness structure, and \code{llrow}, each generated
#' row's log likelihood.
#' @aliases ctStanGenerateFromFit
#' @export
#' @examples
#' gen <- ctGenerateFromFit(ctstantestfit, nsamples=3,fullposterior=TRUE,cores=1)
#' plot(gen$generated$Y[3,,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
ctGenerateFromFit<-function(fit,nsamples=200,fullposterior=FALSE, verboseErrors=FALSE,cores=getOption("mc.cores", 2L),
  intoverstates='fit', effects=c('population','fitted')){
  
  # The julia engine generates from the same forward pass it filters with,
  # drawing each row from its own prior predictive (see R/ctBackendKalman.R and
  # the engine's kalman_trace.jl).
  if(.ctFitIsJulia(fit)){
    return(.ctBackendGenerateFromFit(fit,nsamples=nsamples,
      fullposterior=fullposterior,cores=cores,intoverstates=intoverstates,
      effects=match.arg(effects)))
  }
  if(!'ctStanFit' %in% class(fit)) stop('Not a fit from ctFit().')
  if(!missing(effects)) stop(
    "effects is a backend='julia' option: the generated Stan model draws ",
    'random effects as the fit carried them. Drop the argument, or refit ',
    'with backend="julia".', call.=FALSE)
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
