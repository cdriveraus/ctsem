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

ctFitAutoGetIndividualRestrictions <- function(fits,f0,groupFreeThreshold){
  restrictedPars <- lapply(fits,function(x) x$stanfit$parsteps)
  possibleRestrictions <- unique(unlist(restrictedPars))
  restrictedParsTF <- t(sapply(fits,function(x) possibleRestrictions %in% x$stanfit$parsteps))
  
  restrictedParsMeans <- colMeans(restrictedParsTF)
  ms=f0$setup$matsetup
  names(restrictedParsMeans) <- ms$parname[match(possibleRestrictions,ms$param)]
  
  groupRestrictedPars <- possibleRestrictions[which(restrictedParsMeans > groupFreeThreshold)]
  return(list(groupRestrictedPars=groupRestrictedPars,
    individualProportions=1-restrictedParsMeans, possibleRestrictions=possibleRestrictions))
}


#' @title ctFitAuto
#' @description Fit a ctStan model with automatic parameter selection
#' @param m ctStan model object without time independent predictors. 
#' @param dat Data in long format
#' @param DRIFT Logical, if TRUE, off diagonal drift parameters in the model are tested for inclusion
#' @param DIFFUSION Logical, if TRUE, off diagonal diffusion parameters in the model are tested for inclusion
#' @param fast Logical, if TRUE, do not compute uncertainty hessian / samples in individual level models. 
#' @param initialRestrictions Alternative to the DRIFT / DIFFUSION arguments -- specify explicitly which parameters should be fixed initially, vector of integers based on the $setup$matsetup element of the ctStanFit object, which gives the parameter numbers. Primarily for internal use. 
#' @param individuals Logical, if TRUE, fit individual level models and determine a group model based on the groupFreeThreshold argument.
#' @param groupFreeThreshold Numeric, threshold for group model structure -- if a parameter improves fit in this proportion of individuals or greater, it is freed for all individuals.
#' @param cores Number of CPU cores to use
#' @param ... Additional arguments passed to ctStanFit
#' @details This function is used to automatically select parameters in a ctStan model. Any specified DRIFT / DIFFUSION matrix off diagonals are only included if they significantly improve the likelihood, based on an estimated likelihood ratio test (relying on the Hessian).
#' @return A ctStan fit object
#' @export
#' @importFrom future.apply future_lapply
#' @examples
#' #' \dontrun{
#' testmodel <- ctstantestfit$ctstanmodelbase
#' testmodel$pars$TI1_effect <- NULL
#' testmodel$n.TIpred <- 0
#' testmodel$TIpredNames <- NULL
#' testfit <- ctFitAuto(testmodel, dat = ctstantestdat, DRIFT = TRUE, DIFFUSION = TRUE)
#' summary(testfit)
#' }
ctFitAuto <- function(m, dat, DRIFT=TRUE, DIFFUSION=TRUE,fast=FALSE,initialRestrictions=NA, individuals=FALSE,groupFreeThreshold=.5,cores=2,...){
  
  dat = data.frame(dat)
  
  if(individuals & cores > sum(!duplicated(dat[[m$subjectIDname]]))) cores <- sum(!duplicated(dat[[m$subjectIDname]]))
  if(individuals & !fast) message('Refitting individual level models, this may take some time -- consider specifying fast argument to avoid this!')
  f0 <- ctStanFit(datalong = dat, ctstanmodel = m,fit=F) #fit structure without fitting
  
  ms=f0$setup$matsetup #matrix setup
  ms <- ms[!duplicated(ms$param) & ms$param > 0 & ms$when==0,] #with only estimated parameters
  ms$pass = 1
  if(!is.na(initialRestrictions[1])){ #if initial restrictions are specified, use them
    ms$pass[ms$param %in% initialRestrictions] <- 2
  } else{ 
    if(DRIFT) ms$pass[(ms$matrix==3 & ms$row!=ms$col)] <- 2
    if(DIFFUSION) ms$pass[(ms$matrix==4 & ms$row!=ms$col)] <- 2
    if(!DRIFT && !DIFFUSION) stop("At least one of DRIFT or DIFFUSION must be TRUE")
  }
  
  ffinal=list()
  parsteps = list(ms$param[(ms$pass==2)]) #list containing the parameters to be fixed initially
  #fit to determine which parameters should be free / restricted
  f  <- ctStanFit(datalong = dat, ctstanmodel = f0$ctstanmodelbase, cores=cores,
    optimcontrol=list(parsteps=parsteps,
      carefulfit=ifelse(individuals,FALSE,TRUE),
      parstepsAutoModel=ifelse(individuals,'group',TRUE),
      groupFreeThreshold=groupFreeThreshold,estonly=T),...)
  
  
  parsteps = f$stanfit$parsteps #group / full restrictions
  if(individuals){
    ffinal <- list()
    ffinal$subjectRestrictions <- !f$stanfit$optimfit$subjFreed
    ms=f0$setup$matsetup #setup colnames for subject restriction dataframe output
    colnames(ffinal$subjectRestrictions) <- ms$parname[match(parsteps,ms$param)]
    
    subjparsteps = lapply(1:nrow(f$stanfit$optimfit$subjFreed),function(si){ #get the subject restriction integers
      parsteps[which(f$stanfit$optimfit$subjFreed[si,])]
    })

    ffinal$subjectmodels <- lapply(subjparsteps, function(subjiparsteps){ # get the model structure for each subject
      ctFitAutoRestrictPars(m, f, subjiparsteps)
    })

    if(!fast){
      message('Refitting individual level models')
      if(!requireNamespace("future.apply")) stop("Package \"future.apply\" is required but not installed. Please install it using: install.packages('future.apply')")
        # Save the current future plan and restore it when the function exits
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        future::plan(multisession, workers = cores)
        
        ffinal$subjectfits <- lapply(1:nrow(f$stanfit$optimfit$subjFreed), function(si){
          subjinit <- f$stanfit$optimfit$subjPars[si,]
          if(length(subjparsteps[[si]]) > 0) subjinit <- subjinit[-subjparsteps[[si]] ] #if there are any fixed parameters, remove them from the initial values
          sfit <- ctStanFit(datalong = dat[dat[[m$subjectIDname]]==unique(dat[[m$subjectIDname]])[si],], ctstanmodel = ffinal$subjectmodels[[si]],init=subjinit,optimcontrol=list(carefulfit=F),cores=1,...)
        })#,future.seed=TRUE)
    } #end refit if not fast
  }#end individuals
  
  ffinal$model <- ctFitAutoRestrictPars(m, f, parsteps=f$stanfit$parsteps) #final model structure
  
  #final estimate using restructured model
  
  if(!individuals){
    ffinal <- ctStanFit(datalong = dat, ctstanmodel = ffinal$model,optimcontrol=list(estonly=fast),...)
  ffinal$stanfit$parsteps <- f$stanfit$parsteps #append the parsteps to the final fit
  }
  return(ffinal)
}



#' @title ctFitAutoGroupModel
#' @description Fit a ctStan model with automatic parameter selection for multiple subjects
#' @param m ctStan model object without time independent predictors.
#' @param dat Data in long format
#' @param cores Number of CPU cores to use
#' @param DRIFT Logical, if TRUE, off diagonal drift parameters in the model are tested for inclusion
#' @param DIFFUSION Logical, if TRUE, off diagonal diffusion parameters in the model are tested for inclusion
#' @param groupFreeThreshold Numeric, threshold for group free parameter selection. Default is .5
#' @param ... Additional arguments passed to ctStanFit
#' @details This function is used to automatically select parameters in a ctStan model. Any specified DRIFT / DIFFUSION matrix off diagonals are only included if they significantly improve the likelihood, based on an estimated likelihood ratio test (relying on the Hessian). Subjects are fit one by one, and a group model is determined based on the groupFreeThreshold parameter -- when the proportion of subjects with a parameter free is above this threshold, the parameter is freed in the group model.
#' @return A list containing a list of ctStan fit objects for each subject, and a group model
#' @export
#' @importFrom future.apply future_lapply
#' @examples
#' \dontrun{
#' testmodel <- ctstantestfit$ctstanmodelbase
#' testmodel$pars$TI1_effect <- NULL
#' testmodel$n.TIpred <- 0
#' testmodel$TIpredNames <- NULL
#' testfit <- ctFitAutoGroupModel(testmodel, dat = ctstantestdat, cores=2, DRIFT = TRUE, DIFFUSION = TRUE)
#' ctModelLatex(testfit$groupModel)
#' lapply(testfit$fits,function(x) print(ctStanContinuousPars(x)$DRIFT))
#' }
ctFitAutoGroupModel <- function(m, dat, cores, DRIFT=TRUE, DIFFUSION=TRUE,groupFreeThreshold=.5,...){
  
  message('Fitting individual level models')
  if(!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package \"future.apply\" is required but not installed. Please install it using: install.packages('future.apply')")
  }
  # Save the current future plan and restore it when the function exits
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  future::plan(multisession, workers = cores)
  
  dat=data.frame(dat)
  f0 <- ctStanFit(
    datalong = dat[ dat[[m$subjectIDname]] == dat[[m$subjectIDname]][1],], 
    ctstanmodel = m,fit=F) #fit structure without fitting
  
  fits <- future.apply::future_lapply(unique(dat[[m$subjectIDname]]), function(si){
    datsi <- dat[dat[[m$subjectIDname]]==si,]
    suppressMessages(ctFitAuto(m = m, dat=datsi,DRIFT=DRIFT, DIFFUSION=DIFFUSION, cores=1,fast=TRUE,...))
  },future.seed=TRUE)
  
  groupModelChanges=TRUE
  while(groupModelChanges){
    # Get the parameter restrictions for each subject
    individualRestrictions=ctFitAutoGetIndividualRestrictions(fits,f0,groupFreeThreshold=groupFreeThreshold)
    
    groupRestrictedPars = individualRestrictions$groupRestrictedPars
    
    groupm <- ctFitAutoRestrictPars(m, f0, groupRestrictedPars)
    message('Group model determined, with ', length(groupRestrictedPars), ' candidate parameters fixed and ', 
      length(individualRestrictions$possibleRestrictions) - length(groupRestrictedPars), ' parameters free. Refitting individual level models.')
    fits <- future.apply::future_lapply(unique(dat[[m$subjectIDname]]), function(si){
      datsi <- dat[dat[[m$subjectIDname]]==si,]
      suppressMessages(ctFitAuto(m = m, dat=datsi,initialRestrictions=groupRestrictedPars, cores=1,fast=TRUE,...))
    },future.seed=TRUE)
    
    individualProportions=ctFitAutoGetIndividualRestrictions(fits, f0,groupFreeThreshold=groupFreeThreshold)$individualProportions
    
    if(!any(individualProportions > groupFreeThreshold)) groupModelChanges=FALSE
  }
  
  return(list(fits=fits,groupModel=groupm,individualProportions=individualProportions))
}

