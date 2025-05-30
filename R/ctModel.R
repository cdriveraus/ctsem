#' Define a ctsem model
#' 
#' This function is used to specify a continuous time structural equation model, 
#' which can then be fit to data with function \code{\link{ctStanFit}}.
#' 
#' @param type character string. If 'omx' (default) configures model for maximum likelihood fitting with ctFit, using OpenMx. 
#' If 'ct' or 'dt' configures either continuous ('ct') or discrete ('dt') time 
#' model for Bayesian fitting with \code{\link{ctStanFit}}, using Stan.
#' 
#' @param n.manifest Number of manifest indicators per individual at each measurement occasion / time point.  
#' Manifest variables are included as the first element of the wide data matrix, with all the 1:n.manifest manifest variables 
#' at time 1 followed by those of time 2, and so on.
#' 
#' @param n.latent Number of latent processes.
#' 
#' @param LAMBDA n.manifest*n.latent loading matrix relating latent to manifest variables, 
#' with latent processes 1:n.latent along the columns, and manifest variables
#' 1:n.manifest in the rows.
#' 
#' @param manifestNames n.manifest length vector of manifest variable names as they appear in the data structure, 
#' without any _Tx time point suffix that may be present in wide data.  Defaults to Y1, Y2, etc.
#' 
#' @param manifesttype n.manifest length vector of manifest variable types,defaults to 0 for continuous vars, 1 for binary vars is also possible. 
#' 
#' @param latentNames n.latent length vector of latent variable names 
#' (used for naming parameters, defaults to eta1, eta2, etc).
#' 
#' @param id character string denoting column name containing subject identification variables. 
#' id data may be of any form, though will be coerced internally to an integer sequence rising from 1.
#' @param time character string denoting column name containing timing data. Timing data must be numeric.
#' 
#' @param T0VAR lower triangular n.latent*n.latent cholesky matrix of latent process initial variance / covariance. 
#' "auto" freely estimates all parameters.
#' 
#' @param T0MEANS n.latent*1 matrix of latent process means at first time point, T0. 
#' "auto" freely estimates all parameters.
#' 
#' @param MANIFESTMEANS n.manifest*1 matrix of manifest intercept parameters.  
#' "auto" frees all parameters.
#' 
#' @param MANIFESTVAR lower triangular n.manifest*n.manifest cholesky matrix of variance / covariance 
#' between manifests at each measurement occasion (i.e. measurement error / residual).  
#' "auto" freely estimates variance parameters, 
#' and fixes covariances between manifests to 0. "free" frees all values, including covariances.
#' 
#' @param DRIFT n.latent*n.latent DRIFT matrix of continuous auto and cross effects, 
#' relating the processes over time. 
#' "auto" freely estimates all parameters.
#' 
#' @param CINT n.latent * 1 matrix of latent process intercepts, allowing for non 0 
#' asymptotic levels of the latent processes. Generally only necessary for additional trends and more complex dynamics.
#' 
#' @param DIFFUSION lower triangular n.latent*n.latent cholesky matrix of diffusion process 
#' variance and covariance (latent error / dynamic innovation).
#' "auto" freely estimates all parameters.
#' 
#' @param n.TIpred Number of time independent predictors. 
#' Each TIpredictor is inserted at the right of the data matrix, after the time intervals.
#' 
#' @param TIpredNames n.TIpred length vector of time independent predictor variable names,
#' as they appear in the data structure.  Default names are TI1, TI2, etc.
#' 
#' @param tipredDefault Logical. TRUE sets any parameters with unspecified time independent 
#' predictor effects to have effects estimated, FALSE fixes the effect to zero unless individually specified.
#' 
#' @param n.TDpred Number of time dependent predictor variables in the dataset.  
#' 
#' @param TDpredNames n.TDpred length vector of time dependent predictor variable names, 
#' as they appear in the data structure, without any _Tx time point suffix that may appear in wide data.  
#' Default names are TD1, TD2, etc.
#' 
#' @param Tpoints For type='omx' only. Number of time points, or measurement occasions, in the data.  This will generally be the maximum 
#' number of time points for a single individual, but may be one extra if sample relative time intervals are used, 
#' see \code{\link{ctIntervalise}}. 
#' 
#' @param TRAITVAR For type='omx' only. Either NULL, if no trait / unobserved heterogeneity effect, 
#' or lower triangular n.latent*n.latent cholesky matrix of trait variance / covariance across subjects.
#' "auto" freely estimates all parameters.
#' 
#' @param T0TRAITEFFECT For type='omx' only. Either NULL, if no trait / individual heterogeneity effect, 
#' or lower triangular n.latent*n.latent cholesky matrix of initial trait variance / covariance.
#' "auto" freely estimates all parametrers, if the TRAITVAR matrix is specified.
#' 
#' @param MANIFESTTRAITVAR For type='omx' only. Either NULL (default) if no trait variance / individual heterogeneity in the level of
#' the manifest indicators, otherwise a lower triangular n.manifest * n.manifest variance / covariance matrix. 
#' Set to "auto" to include and free all parameters - but identification problems will arise if \code{TRAITVAR} is 
#' also set.
#' 
#' @param TDPREDMEANS For type='omx' only. (n.TDpred * (Tpoints - 1)) rows * 1 column matrix of time dependent predictor means.
#' If 'auto', the means are freely estimated.  Otherwise, 
#' the means for the Tpoints observations of your first time dependent predictor 
#' are followed by those of TDpred 2, and so on.
#' 
#' @param TDPREDEFFECT n.latent*n.TDpred matrix of effects from time dependent predictors to latent processes.
#' Effects from 1:n.TDpred columns TDpredictors go to 1:n.latent rows of latent processes.
#' "auto" freely estimates all parameters.
#' 
#' @param T0TDPREDCOV For type='omx' only. n.latent rows * (Tpoints * n.TDpred) columns covariance matrix 
#' between latents at T0 and time dependent predictors.
#' Default of "auto" restricts covariance to 0, which is consistent with covariance to other time points. 
#' To freely estimate parameters, specify either 'free', or the desired matrix.
#' 
#' @param TDPREDVAR For type='omx' only. lower triangular (n.TDpred * Tpoints) rows 
#' * (n.TDpred * Tpoints) columns variance / covariance
#' cholesky matrix for time dependent predictors.
#' "auto" (default) freely estimates all parameters.
#' 
#' @param TRAITTDPREDCOV For type='omx' only. n.latent rows * (n.TDpred*Tpoints) columns covariance matrix of 
#' latent traits and time dependent predictors. Defaults to zeroes, 
#' assuming predictors are independent of subjects baseline levels. When predictors depend on the subjects,
#' this should instead be set to 'free' or manually specified.
#' The Tpoints columns of the first preditor are followed by those of the second and so on.
#' Covariances with the trait variance of latent process 1 are specified in row 1, process 2 in row 2, etc.
#' "auto" (default) sets this matrix to zeroes, (if both traits and time dependent predictors exist, otherwise this matrix is set to NULL, and ignored in any case).
#' 
#' @param TDTIPREDCOV For type='omx' only. (n.TDpred * Tpoints) rows * n.TIpred columns covariance
#' matrix between time dependent and time independent predictors.
#' "auto" (default) freely estimates all parameters.
#'  
#' @param TIPREDMEANS For type='omx' only. n.TIpred * 1 matrix of time independent predictor means.
#' If 'auto', the means are freely estimated.  
#' 
#' @param TIPREDEFFECT For type='omx' only. n.latent*n.TIpred effect matrix of time independent predictors on latent processes.
#' "auto" freely estimates all parameters and generates starting values. TIPREDEFFECT parameters for type='stan' are estimated
#' by default on all subject level parameters, to restrict this, 
#' manually edit the model object after creation.
#' 
#' @param T0TIPREDEFFECT For type='omx' only.n.latent*n.TIpred effect matrix of time independent 
#' predictors on latents at T0. "auto" freely estimates all parameters, though note that under the default 
#' setting of \code{stationary} for ctFit, this matrix is ignored as the effects are determined based on
#' the overall process parameters.
#' 
#' @param TIPREDVAR For type='omx' only.lower triangular n.TIpred * n.TIpred Cholesky decomposed covariance
#' matrix for all time independent predictors.
#' "auto" (default) freely estimates all parameters.
#' 
#' @param PARS for types 'ct' and 'dt' only. May be of any structure, only needed to contain extra parameters for certain non-linear models.
#' 
#' @param startValues For type='omx' only. A named vector, where the names of each value must match a parameter in the specified model,
#' and the value sets the starting value for that parameter during optimization.
#' If not set, random starting values representing relatively stable processes with small effects and 
#' covariances are generated by ctFit.  
#' Better starting values may improve model fit speed and the chance of an appropriate model fit. 
#' 
#' @param silent Suppress all output to console. 
#' 
#' @examples 
#'  ### Frequentist example:
#'  ### impulse and level change time dependent predictor 
#'  ### example from Driver, Oud, Voelkle (2015)
#'  data('ctExample2')
#'  tdpredmodel <- ctModel(n.manifest = 2, n.latent = 3, n.TDpred = 1, 
#'  Tpoints = 8, manifestNames = c('LeisureTime', 'Happiness'), 
#'  TDpredNames = 'MoneyInt', 
#'  latentNames = c('LeisureTime', 'Happiness', 'MoneyIntLatent'),
#'  LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3), TRAITVAR = "auto")
#'
#'  tdpredmodel$TRAITVAR[3, ] <- 0
#'  tdpredmodel$TRAITVAR[, 3] <- 0
#'  tdpredmodel$DIFFUSION[, 3] <- 0
#'  tdpredmodel$DIFFUSION[3, ] <- 0
#'  tdpredmodel$T0VAR[3, ] <- 0
#'  tdpredmodel$T0VAR[, 3] <- 0
#'  tdpredmodel$CINT[3] <- 0
#'  tdpredmodel$T0MEANS[3] <- 0
#'  tdpredmodel$TDPREDEFFECT[3, ] <- 1
#'  tdpredmodel$DRIFT[3, ] <- 0
#'  
#'  
#' ###Bayesian example:
#' model<-ctModel(type='ct',
#' n.latent=2, latentNames=c('eta1','eta2'),
#' n.manifest=2, manifestNames=c('Y1','Y2'),
#' n.TDpred=1, TDpredNames='TD1', 
#' n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
#' LAMBDA=diag(2))
#'
#' 
#' @export

ctModel<-function(LAMBDA, type='omx',n.manifest = 'auto', n.latent='auto', Tpoints=NULL, 
  manifestNames='auto', manifesttype=rep(0,nrow(LAMBDA)),latentNames='auto', id='id',time='time', silent=FALSE,
  T0VAR="auto", T0MEANS="auto", MANIFESTMEANS="auto", MANIFESTVAR="diag", 
  DRIFT="auto", CINT=0, DIFFUSION="auto",
  n.TDpred='auto', TDpredNames='auto', 
  n.TIpred='auto', TIpredNames='auto', tipredDefault=TRUE,
  TRAITVAR=NULL, T0TRAITEFFECT=NULL,
  MANIFESTTRAITVAR=NULL,
  TDPREDMEANS="auto", TDPREDEFFECT="auto", 
  T0TDPREDCOV="auto", TDPREDVAR="auto",  
  TRAITTDPREDCOV="auto", TDTIPREDCOV='auto',
  TIPREDMEANS="auto", TIPREDEFFECT="auto", 
  T0TIPREDEFFECT="auto", TIPREDVAR="auto", 
  PARS=NULL,
  startValues=NULL){
  
  #get dimensions
  if(is.null(n.manifest) || is.null(n.latent) || all(n.manifest %in% 'auto') || all(n.latent %in% 'auto')){
    if(is.matrix(LAMBDA)){
      if(!silent) message('System dimensions inferred from LAMBDA')
      n.manifest <- nrow(LAMBDA)
      n.latent <- ncol(LAMBDA)
    } else {
      if(!all(manifestNames %in% 'auto') && !all(latentNames %in% 'auto')){
        if(!silent) message('System dimensions inferred from manifestNames and latentNames')
        n.manifest <- length(manifestNames)
        n.latent <- length(latentNames)
      } else stop('LAMBDA must either a matrix, n.manifest and n.latent must be specified, or manifestNames and latentNames must be specified!')
    }
  }
  
  if(is.null(n.TDpred) || all(n.TDpred %in% 'auto')){
    if(is.matrix(TDPREDEFFECT)){
      if(!silent) message('n.TDpred inferred from TDPREDEFFECT')
      n.TDpred <- ncol(TDPREDEFFECT)
    } else {
      if(!all(TDpredNames %in% 'auto')){
        if(!silent) message('n.TDpred inferred  inferred from TDpredNames')
        n.TDpred <- length(TDpredNames)
      } else n.TDpred <- 0
    } 
  }
  
  if(is.null(n.TIpred) || all(n.TIpred %in% 'auto')){
    if(!all(TIpredNames %in% 'auto')){
      if(!silent) message('n.TIpred inferred  inferred from TIpredNames')
      n.TIpred <- length(TIpredNames)
    } else n.TIpred <- 0
  } 
  
  
  
  #####FUNCTIONS
  
  
  checkSymmetry<-function(x){  #this checks the symmetry of matrix x, if not symmetric it stops and outputs an error
    if(isSymmetric(x)==F){
      stop(paste0(substitute(x)," must be symmetric"))
    }
  }
  
  
  if(type=='omx') message('Type "omx" is still supported but requires ctsemOMX package installation. "ct" or "dt" are recommended types.')
  if(type=='omx' & is.null(Tpoints)) stop('Type "omx" requires Tpoints specified!')
  
  
  
  
  ###### RUN SEQUENCE
  
  if(type!='omx' && is.null(Tpoints)) Tpoints<-3
  # if(type=='omx' && is.null(Tpoints)) stop('Tpoints must be specified for type="omx"')

  
  #names
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  if(all(latentNames=='auto')) latentNames=paste0('eta',1:n.latent)
  
  if(length(manifestNames) != n.manifest) stop(paste0(
    "Length of manifestNames (",length(manifestNames),") does not equal n.manifest (",n.manifest,")!")) 
  if(length(latentNames) != n.latent) stop(paste0(
    "Length of latentNames (",length(latentNames),") does not equal n.latent(",n.latent,")!"))
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop(paste0(
      "Length of TDpredNames (",length(TDpredNames),") does not equal n.TDpred(",n.TDpred,")!")) 
  } else {
    TDpredNames=c()
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop(paste0(
      "Length of TIpredNames (",length(TIpredNames),') does not equal n.TIpred(',n.TIpred,')!')) 
  } else {
    TIpredNames=c()
  }
  
  mats <- ctStanMatricesList()
  for(m in names(mats$base)){
    if(!is.null(get(m))){ #if the matrix is specified
      val <- get(m)
      if(!all(val %in% c('auto','diag')) && !is.matrix(val)){ #and not auto or a matrix
        if(m %in% 'PARS') mat <- matrix(PARS,ncol=1) else {
          mat<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames, #set the basic structure
            manifestNames=manifestNames,latentNames=latentNames,matrixname=m,n.latent=n.latent,
            n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
        }
        counter <- 0 #then start filling with vector
        if(length(val) != length(mat) && length(val) !=1) stop(paste0(m,' needs ',length(mat),' values for ',nrow(mat),' * ',ncol(mat),' matrix, but has ',length(val)))
        if(length(val)==1 && length(mat) > 1) {
          if(!silent) message(m,' specified via single value -- filling ',nrow(mat),' * ',ncol(mat),' matrix:')
          val <- rep(val,length(mat))
        } else if(!silent) message(paste0(m,' vector spec input rowwise into ',nrow(mat),' * ',ncol(mat),' matrix:'))
        for(ri in 1:nrow(mat)){
          for(ci in 1:ncol(mat)){
            counter <- counter + 1
            mat[ri,ci] <- val[counter]
          }
        }
        
        if(!silent) print(mat,right=TRUE)
        assign(m, mat)
      }
    }
  }
  
  
  #matrices
  if(T0MEANS[1]=="auto") T0MEANS<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="T0MEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(ncol(T0MEANS)>1) stop("Specified T0MEANS matrix has more than one column")
  if(nrow(T0MEANS)!=n.latent) stop ("Specified T0MEANS matrix rows not equal to n.latent")
  rownames(T0MEANS)=latentNames
  
  
  
  
  
  if(T0VAR[1]=="auto") T0VAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="T0VAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(nrow(T0VAR)!=n.latent | ncol(T0VAR)!=n.latent) stop("Dimensions of T0VAR matrix are not n.latent * n.latent")
  dimnames(T0VAR)=list(latentNames,latentNames)
  
  
  
  
  
  if(MANIFESTMEANS[1]=="auto") MANIFESTMEANS<- ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTMEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(ncol(MANIFESTMEANS)>1) stop("Specified MANIFESTMEANS matrix has more than one column")
  if(nrow(MANIFESTMEANS)!=n.manifest) stop ("Specified MANIFESTMEANS matrix rows not equal to n.manifest")
  rownames(MANIFESTMEANS)=manifestNames
  
  
  

    if(all(MANIFESTVAR %in% "auto")) {
      MANIFESTVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
        manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTVAR",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    }
    
    if(all(MANIFESTVAR %in% 'diag')){
      MANIFESTVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
        manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTVAR",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
      MANIFESTVAR[row(MANIFESTVAR) != col(MANIFESTVAR)] <- 0
    }
    
    if(nrow(MANIFESTVAR)!=n.manifest | ncol(MANIFESTVAR)!=n.manifest) stop("Dimensions of MANIFESTVAR matrix are not n.manifest * n.manifest")
  
  dimnames(MANIFESTVAR)=list(manifestNames,manifestNames)
  
  
  
  if(all(DRIFT %in% "auto")) DRIFT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DRIFT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(all(DRIFT %in% 'diag')) {
    DRIFT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DRIFT",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    DRIFT[row(DRIFT) != col(DRIFT)] <- 0
  }
  if(nrow(DRIFT)!= n.latent | ncol(DRIFT)!= n.latent) stop("Dimensions of DRIFT matrix are not n.latent * n.latent")
  dimnames(DRIFT)=list(latentNames,latentNames)
  
  
  
  
  if(CINT[1]=="auto") CINT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="CINT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(ncol(CINT)!=1 | nrow(CINT)!=n.latent) stop("Dimensions of CINT are not n.latent * 1")
  rownames(CINT)=latentNames
  
  
  
  
  if(all(DIFFUSION %in% "auto")) DIFFUSION<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DIFFUSION",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(all(DIFFUSION %in% 'diag')) {
    DIFFUSION<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DIFFUSION",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    DIFFUSION[row(DIFFUSION) != col(DIFFUSION)] <- 0
  }
  if(nrow(DIFFUSION)!=n.latent | ncol(DIFFUSION) != n.latent) stop("Dimensions of DIFFUSION are not n.latent * n.latent")
  dimnames(DIFFUSION)=list(latentNames,latentNames)
  
  
  
  
  
  if(!is.null(TRAITVAR)){
    
    if(all(TRAITVAR=="auto")) TRAITVAR<-ctLabel(TDpredNames=TDpredNames,
      TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TRAITVAR",
      n.latent=n.latent,n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    
    if(nrow(TRAITVAR)+ncol(TRAITVAR)!=n.latent*2) stop ("
    Dimensions of TRAITVAR are not n.latent * n.latent")
    dimnames(TRAITVAR)=list(latentNames,latentNames)
    
    if(is.null(T0TRAITEFFECT) || all(T0TRAITEFFECT == 'auto')) T0TRAITEFFECT<-
      ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,
        latentNames=latentNames,matrixname="T0TRAITEFFECT",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    dimnames(T0TRAITEFFECT)=list(latentNames,latentNames)
    
    if(!is.null(T0TRAITEFFECT) && nrow(T0TRAITEFFECT)+ncol(T0TRAITEFFECT)!=n.latent*2) stop ("Dimensions of TRAITVAR are not n.latent * n.latent")
    
  }
  
  
  
  
  
  
  
  
  if(!is.null(MANIFESTTRAITVAR)){
    
    if(MANIFESTTRAITVAR[1]=="auto") MANIFESTTRAITVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTTRAITVAR",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    
    if(nrow(MANIFESTTRAITVAR)+ncol(MANIFESTTRAITVAR)!=n.manifest*2) stop ("Dimensions of MANIFESTTRAITVAR are not n.manifest * n.manifest")
    #   if(is.null(MANIFESTTRAITVAR)) MANIFESTTRAITVAR <- matrix(0,nrow=n.manifest,ncol=n.manifest)
    dimnames(MANIFESTTRAITVAR)=list(manifestNames,manifestNames)
  }
  
  
  
  
  #effect matrix of time dependent predictors (with time independent effect strength)
  if(n.TDpred > 0){
    if(all(TDPREDEFFECT=="auto"))  {
      
      TDPREDEFFECT<-ctLabel(TDpredNames=TDpredNames,
        TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDEFFECT",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
      
      if(nrow(TDPREDEFFECT) != n.latent | ncol(TDPREDEFFECT) != n.TDpred) stop("Dimensions of TDPREDEFFECT are not n.latent * n.TDpred")
    }
    dimnames(TDPREDEFFECT)=list(latentNames,TDpredNames)
    
    # means of TD predictors
    if(all(TDPREDMEANS=="auto") & n.TDpred > 0) TDPREDMEANS <- ctLabel(TDpredNames=TDpredNames,
      TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDMEANS",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    
    if(nrow(TDPREDMEANS) != n.TDpred*(Tpoints) | ncol(TDPREDMEANS) != 1) stop(
      "Dimensions of TDPREDMEANS are not (n.TDpred*Tpoints) rows * 1 column")
  }
  
  if(n.TDpred==0) {
    TDPREDEFFECT<-NULL
    TDPREDMEANS<-NULL
  }
  
  
  # means of TI predictors
  if(TIPREDMEANS[1]=="auto" & n.TIpred >0) TIPREDMEANS <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TIPREDMEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(n.TIpred==0) TIPREDMEANS<-NULL
  
  if(n.TIpred > 0) if(nrow(TIPREDMEANS) != n.TIpred | ncol(TIPREDMEANS) != 1) stop(
    "Dimensions of TIPREDMEANS are not n.TIpred * 1")
  
  
  
  
  if(all(T0TDPREDCOV=="free") && n.TDpred>0) T0TDPREDCOV<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="T0TDPREDCOV",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(all(T0TDPREDCOV=="auto") && n.TDpred>0) T0TDPREDCOV<-matrix(0,n.latent,n.TDpred*(Tpoints))
  if(n.TDpred==0) T0TDPREDCOV<-NULL
  
  if(n.TDpred>0) if(nrow(T0TDPREDCOV) != n.latent | ncol(T0TDPREDCOV) != n.TDpred*(Tpoints))  stop(
    "Dimensions of T0TDPREDCOV are not n.latent * (n.TDpred*Tpoints)")
  
  
  
  
  
  if(all(TIPREDEFFECT=="auto") && n.TIpred>0) TIPREDEFFECT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TIPREDEFFECT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(n.TIpred==0) TIPREDEFFECT<-NULL
  if(n.TIpred>0) if(nrow(TIPREDEFFECT) != n.latent | ncol(TIPREDEFFECT) != n.TIpred) stop("Dimensions of TIPREDEFFECT are not n.latent * n.TIpred")
  
  
  
  
  
  if(all(T0TIPREDEFFECT=="auto") && n.TIpred > 0) T0TIPREDEFFECT <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="T0TIPREDEFFECT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(n.TIpred==0) T0TIPREDEFFECT <- NULL
  
  if(n.TIpred>0) if(nrow(T0TIPREDEFFECT) != n.latent | ncol(T0TIPREDEFFECT) != n.TIpred) stop("Dimensions of T0TIPREDEFFECT are not n.latent * n.TIpred")
  
  
  
  
  
  if(all(TIPREDVAR=="auto") && n.TIpred > 0)  TIPREDVAR <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TIPREDVAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints) 
  
  if(n.TIpred==0) TIPREDVAR <- NULL
  
  if(n.TIpred > 0) if(nrow(TIPREDVAR) != n.TIpred | ncol(TIPREDVAR) != n.TIpred) stop(
    "Dimensions of TIPREDVAR are not n.TIpred * n.TIpred")
  
  
  
  
  
  if(all(TDPREDVAR=="auto") && n.TDpred >0)  TDPREDVAR <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDVAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)  
  
  if(n.TDpred==0) TDPREDVAR<-NULL
  
  if(n.TDpred > 0) if(nrow(TDPREDVAR) != n.TDpred*(Tpoints) | ncol(TDPREDVAR) != n.TDpred*(Tpoints)) stop(
    "Dimensions of TDPREDVAR are not (n.TDpred*Tpoints) * (n.TDpred*Tpoints)")
  
  
  
  
  if(all(TDTIPREDCOV=="free") & n.TDpred > 0 & n.TIpred > 0)  TDTIPREDCOV <- ctLabel(TDpredNames=TDpredNames,
    TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDTIPREDCOV",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(all(TDTIPREDCOV=="auto") & n.TDpred > 0 & n.TIpred > 0)  TDTIPREDCOV <- matrix(0,(n.TDpred*(Tpoints)), n.TIpred)
  if(n.TDpred == 0 | n.TIpred == 0) TDTIPREDCOV<-NULL
  if(n.TDpred > 0 & n.TIpred > 0) if(nrow(TDTIPREDCOV) != (n.TDpred*(Tpoints)) | ncol(TDTIPREDCOV) != n.TIpred) stop(
    "Dimensions of TDTIPREDCOV are not (n.TDpred*Tpoints) * n.TIpred")
  
  
  
  
  
  # TDpredictor and trait covariances
  if(n.TDpred > 0) { #then set covariance matrix for heterogeneity and predictors    
    if(all(TRAITTDPREDCOV == "free")) TRAITTDPREDCOV<-ctLabel(TDpredNames=TDpredNames,
      TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TRAITTDPREDCOV",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)  
    if(all(TRAITTDPREDCOV == "auto")) TRAITTDPREDCOV<-matrix(0,n.latent, (n.TDpred*(Tpoints)))
    if(nrow(TRAITTDPREDCOV) != n.latent | ncol(TRAITTDPREDCOV) != n.TDpred*(Tpoints)) stop(
      "Dimensions of TRAITTDPREDCOV are not n.latent * (n.TDpred*Tpoints)")
  }
  
  if(n.TDpred==0)  TRAITTDPREDCOV<-NULL #if no predictors, assume no covariance between heterogeneity and predictors
  
  
  
  ###model checks
  
  #lower triangular check
  for(tempmatname in c('T0VAR','MANIFESTVAR','TRAITVAR','MANIFESTTRAITVAR','TDPREDVAR','DIFFUSION','TIPREDVAR')){
    assign('tempmat',get(tempmatname))
    if(is.null(tempmat)) next
    if(all(suppressWarnings(as.numeric(tempmat[upper.tri(tempmat)])) %in% 0)) next
    warning(paste0(tempmatname,' is not lower triangular! Covariance type matrices should usually be specified in the appropriate lower-triangular form.'))
  }
  
  # mvaroffdiag=MANIFESTVAR[!diag(1,nrow(MANIFESTVAR))]
  # if(nrow(MANIFESTVAR) > 1 && !all(mvaroffdiag %in% 0)) stop('MANIFESTVAR should be diagonal!')
  
  
  if(any(dim(LAMBDA)!=c(n.manifest,n.latent))) stop("Incorrect LAMBDA structure specified - check number or rows and columns")
  dimnames(LAMBDA)=list(manifestNames,latentNames)
  
  sapply(c(manifestNames,latentNames,TDpredNames,TIpredNames, time,id),function(x){
    if(grepl('\\W',x)) stop(paste0(x,' contains symbols, variable names must be alphanumerics only please!'))
  })
  
  for(names in c('manifestNames','TIpredNames','TDpredNames','latentNames')){
    if(any(duplicated(get(names)))) stop(paste0('Duplicate names in ',names))
  }
  
  if(any(manifesttype>0 ) && all(CINT %in% 0)) warning('CINT usually needs to be specified for non-continuous variables -- consider fixing relevant MANIFESTMEANS to zero instead')
  
  
  
  
  # completemodel<-list(n.manifest,n.latent,n.TDpred,n.TIpred,Tpoints,LAMBDA,
  #   manifestNames,latentNames,TDpredNames,TIpredNames,
  #   T0VAR,T0MEANS,MANIFESTMEANS,MANIFESTVAR,DRIFT,CINT,DIFFUSION,
  #   TRAITVAR,T0TRAITEFFECT,MANIFESTTRAITVAR,
  #   TDPREDEFFECT, TDPREDMEANS, T0TDPREDCOV, TRAITTDPREDCOV,
  #   TIPREDEFFECT, TIPREDMEANS, T0TIPREDEFFECT,
  #   TIPREDVAR, TDPREDVAR,TDTIPREDCOV,  PARS, 
  #   startValues,id,time,manifestType)
  # 
  # names(completemodel)<-c("n.manifest","n.latent","n.TDpred","n.TIpred","Tpoints","LAMBDA",
  #   "manifestNames","latentNames","TDpredNames","TIpredNames",
  #   "T0VAR","T0MEANS","MANIFESTMEANS","MANIFESTVAR","DRIFT","CINT","DIFFUSION",
  #   "TRAITVAR",'T0TRAITEFFECT',"MANIFESTTRAITVAR",
  #   'TDPREDEFFECT', 'TDPREDMEANS', 'T0TDPREDCOV', 'TRAITTDPREDCOV',
  #   'TIPREDEFFECT', 'TIPREDMEANS', 'T0TIPREDEFFECT',
  #   'TIPREDVAR', 'TDPREDVAR','TDTIPREDCOV',  'PARS',
  #   "startValues","id","time","manifestType")
  
  #rewrite last 2 pieces of code into a single element:
  completemodel <- list(`n.manifest`=n.manifest, `n.latent`=n.latent, `n.TDpred`=n.TDpred, `n.TIpred`=n.TIpred, 
    `Tpoints`=Tpoints, `LAMBDA`=LAMBDA, `manifestNames`=manifestNames, `latentNames`=latentNames, 
    `TDpredNames`=TDpredNames, `TIpredNames`=TIpredNames, `T0VAR`=T0VAR, `T0MEANS`=T0MEANS, 
    `MANIFESTMEANS`=MANIFESTMEANS, `MANIFESTVAR`=MANIFESTVAR, `DRIFT`=DRIFT, `CINT`=CINT, 
    `DIFFUSION`=DIFFUSION, `TRAITVAR`=TRAITVAR, `T0TRAITEFFECT`=T0TRAITEFFECT, 
    `MANIFESTTRAITVAR`=MANIFESTTRAITVAR, `TDPREDEFFECT`=TDPREDEFFECT, 
    `TDPREDMEANS`=TDPREDMEANS, `T0TDPREDCOV`=T0TDPREDCOV, `TRAITTDPREDCOV`=TRAITTDPREDCOV, 
    `TIPREDEFFECT`=TIPREDEFFECT, `TIPREDMEANS`=TIPREDMEANS, `T0TIPREDEFFECT`=T0TIPREDEFFECT, 
    `TIPREDVAR`=TIPREDVAR, `TDPREDVAR`=TDPREDVAR, `TDTIPREDCOV`=TDTIPREDCOV, `PARS`=PARS, 
    `startValues`=startValues, `id`=id, `time`=time, `manifesttype`=manifesttype)
  
  
  
  
  if(type=='omx') class(completemodel)<-"ctsemInit"
  
  if(type %in% c('ct','dt','stanct','standt')) completemodel<-ctStanModel(completemodel,type=type,tipredDefault= tipredDefault)
  if(!type %in% c('omx','ct','dt','stanct','standt')) stop('type must be either omx, ct, dt!')
  
  return(completemodel)
}
