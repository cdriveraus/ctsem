#' @title ctFitAuto
#' @description Fit a ctStan model with automatic parameter selection
#' @param m ctStan model object without time independent predictors. 
#' @param dat Data in long format
#' @param DRIFT Logical, if TRUE, off diagonal drift parameters in the model are tested for inclusion
#' @param DIFFUSION Logical, if TRUE, off diagonal diffusion parameters in the model are tested for inclusion
#' @param ... Additional arguments passed to ctStanFit
#' @details This function is used to automatically select parameters in a ctStan model. Any specified DRIFT / DIFFUSION matrix off diagonals are only included if they significantly improve the likelihood, based on an estimated likelihood ratio test (relying on the Hessian).
#' @return A ctStan fit object
#' @export
#' @examples
#' #' \dontrun{
#' testmodel <- ctstantestfit$ctstanmodelbase
#' testmodel$pars$TI1_effect <- NULL
#' testmodel$n.TIpred <- 0
#' testmodel$TIpredNames <- NULL
#' testfit <- ctFitAuto(testmodel, dat = ctstantestdat, DRIFT = TRUE, DIFFUSION = TRUE)
#' summary(testfit)
#' }
ctFitAuto <- function(m, dat, DRIFT=TRUE, DIFFUSION=TRUE,...){
  
  f0 <- ctStanFit(datalong = dat, ctstanmodel = m,fit=F) #fit structure without fitting
  ms=f0$setup$matsetup #matrix setup
  ms <- ms[!duplicated(ms$param) & ms$param > 0 & ms$when==0,] #with only estimated parameters
  ms$pass = 1
  if(DRIFT) ms$pass[(ms$matrix==3 & ms$row!=ms$col)] <- 2
  if(DIFFUSION) ms$pass[(ms$matrix==4 & ms$row!=ms$col)] <- 2
  if(!DRIFT && !DIFFUSION) stop("At least one of DRIFT or DIFFUSION must be TRUE")
  
  parsteps = list(ms$param[(ms$pass==2)]) #list containing the parameters to be fixed initially
  
  #fit to determine which parameters should be free
  f  <- ctStanFit(datalong = dat, ctstanmodel = m,optimcontrol=list(parsteps=parsteps,parstepsAutoModel=T,estonly=T),...)
  
  #final model structure:
  mfinal <- m
  ms=data.frame(f$standata$matsetup)
  
  for(fixedpari in f$stanfit$parsteps){
    par_ms <- ms[ms$param == fixedpari,,drop=FALSE][1,]
    parmat <- names(sort(ctStanMatricesList()$all))[par_ms$matrix]
    parrow <- par_ms$row
    parcol <- par_ms$col
    mfinal$pars[mfinal$pars$matrix== parmat & mfinal$pars$row==parrow & mfinal$pars$col==parcol,'value'] <- 0
  }
  
  #final estimate
  ffinal <- ctStanFit(datalong = dat, ctstanmodel = mfinal,...)
  return(ffinal)
}
