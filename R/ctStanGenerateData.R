#' Add a $generated object to ctstanfit object, with random data generated from posterior of ctstanfit object
#'
#' @param fit ctstanfit object
#'
#' @return Matrix of generated data -- one dataset per iteration, according to original time and missingness structure.
#' @export
#'
#' @examples
#' gen <- ctStanGenerateData(ctstantestfit)
#' plot(gen$generated$Y[3,,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
ctStanGenerateData<-function(fit){
  if(class(fit)!='ctStanFit') stop('Not a ctStanFit object!')
  if(class(fit$stanfit)!='stanfit') {
    umat=t(fit$stanfit$rawposterior)
  } else  {
    umat <- stan_unconstrainsamples(fit$stanfit,fit$standata)
  }
  
  if(fit$setup$recompile) {
    message('Compilation needed -- compiling (roughly 1 min)')
    genm <- rstan::stan_model(model_code = ctStanModelWriter(ctstanmodel = fit$ctstanmodel,gendata = TRUE,extratforms = fit$setup$extratforms))
  } else {
    genm <- stanmodels$ctsmgen
  }
  message('Generating data from posterior')
  genf <- stan_reinitsf(genm,fit$standata) 
  fit$generated$Y <- array(apply(umat, 2, function(x) rstan::constrain_pars(genf,x)$Y),dim=c(nrow(fit$standata$Y),fit$ctstanmodel$n.manifest,ncol(umat)))
  fit$generated$Y <- aperm(fit$generated$Y,c(1,3,2))
  fit$generated$Y[fit$generated$Y==99999] <- NA
  return(fit)
}
