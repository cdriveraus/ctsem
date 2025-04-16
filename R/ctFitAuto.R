ctFitAutoRestrictPars <- function(m, f, parsteps){
  ms=data.frame(f$standata$matsetup)
  for(fixedpari in parsteps){
    par_ms <- ms[ms$param == fixedpari,,drop=FALSE][1,]
    parmat <- names(sort(ctStanMatricesList()$all))[par_ms$matrix]
    parrow <- par_ms$row
    parcol <- par_ms$col
    m$pars[m$pars$matrix== parmat & m$pars$row==parrow & m$pars$col==parcol,'value'] <- 0
    m$pars[m$pars$matrix== parmat & m$pars$row==parrow & m$pars$col==parcol,'param'] <- NA
    m$pars[m$pars$matrix== parmat & m$pars$row==parrow & m$pars$col==parcol,'transform'] <- NA
    m$pars[m$pars$matrix== parmat & m$pars$row==parrow & m$pars$col==parcol,'indvarying'] <- FALSE
  }
  return(m)
}


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
  #fit to determine which parameters should be free / restricted
  f  <- ctStanFit(datalong = dat, ctstanmodel = m,optimcontrol=list(parsteps=parsteps,parstepsAutoModel=T,estonly=T),...)
  
  mfinal <- ctFitAutoRestrictPars(m, f, parsteps=f$stanfit$parsteps) #final model structure
  
  #final estimate using restructured model
  ffinal <- ctStanFit(datalong = dat, ctstanmodel = mfinal,...)
  ffinal$stanfit$parsteps <- f$stanfit$parsteps #append the parsteps to the final fit
  return(ffinal)
}



#' @title ctFitAutoGroupModel
#' @description Fit a ctStan model with automatic parameter selection for multiple subjects
#' @param m ctStan model object without time independent predictors.
#' @param dat Data in long format
#' @param DRIFT Logical, if TRUE, off diagonal drift parameters in the model are tested for inclusion
#' @param DIFFUSION Logical, if TRUE, off diagonal diffusion parameters in the model are tested for inclusion
#' @param groupFreeThreshold Numeric, threshold for group free parameter selection. Default is .5
#' @param ... Additional arguments passed to ctStanFit
#' @details This function is used to automatically select parameters in a ctStan model. Any specified DRIFT / DIFFUSION matrix off diagonals are only included if they significantly improve the likelihood, based on an estimated likelihood ratio test (relying on the Hessian). Subjects are fit one by one, and a group model is determined based on the groupFreeThreshold parameter -- when the proportion of subjects with a parameter free is above this threshold, the parameter is freed in the group model.
#' @return A list containing a list of ctStan fit objects for each subject, and a group model
#' @export
#' @examples
#' \dontrun{
#' testmodel <- ctstantestfit$ctstanmodelbase
#' testmodel$pars$TI1_effect <- NULL
#' testmodel$n.TIpred <- 0
#' testmodel$TIpredNames <- NULL
#' testfit <- ctFitAutoGroupModel(testmodel, dat = ctstantestdat, DRIFT = TRUE, DIFFUSION = TRUE)
#' ctModelLatex(testfit$groupModel)
#' lapply(testfit$fits,function(x) print(ctStanContinuousPars(x)$DRIFT))
#' }
ctFitAutoGroupModel <- function(m, dat, DRIFT=TRUE, DIFFUSION=TRUE,groupFreeThreshold=.5, ...){
  
  
  dat=data.frame(dat)
  f0 <- ctStanFit(
    datalong = dat[ dat[[m$subjectIDname]] == dat[[m$subjectIDname]][1],], 
    ctstanmodel = m,fit=F) #fit structure without fitting
  
  message('Fitting individual level models')
  fits <- future.apply::future_lapply(unique(dat[[m$subjectIDname]]), function(si){
    datsi <- dat[dat[[m$subjectIDname]]==si,]
    ctFitAuto(m = m, dat=datsi,DRIFT=DRIFT, DIFFUSION=DIFFUSION, cores=1)
  },future.seed=TRUE)

  restrictedPars <- lapply(fits,function(x) x$stanfit$parsteps)
  possibleRestrictions <- unique(unlist(restrictedPars))
  restrictedParsTF <- t(sapply(fits,function(x) possibleRestrictions %in% x$stanfit$parsteps))
  
  restrictedParsMeans <- colMeans(restrictedParsTF)
  
  groupRestrictedPars <- possibleRestrictions[which(restrictedParsMeans > groupFreeThreshold)]

  groupm <- ctFitAutoRestrictPars(m, f0, groupRestrictedPars)
  message('Group model determined, with ', length(groupRestrictedPars), ' parameters fixed and ', 
          length(possibleRestrictions) - length(groupRestrictedPars), ' parameters free. Refitting individual level models.')
  refits <- future.apply::future_lapply(unique(dat[[m$subjectIDname]]), function(si){
    datsi <- dat[dat[[m$subjectIDname]]==si,]
    ctFitAuto(m = groupm, dat=datsi,DRIFT=DRIFT, DIFFUSION=DIFFUSION, cores=1)
  },future.seed=TRUE)
  
  return(list(fits=refits,groupModel=groupm))
}
