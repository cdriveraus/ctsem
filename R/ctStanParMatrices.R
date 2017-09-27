#' Returns model equation matrices from a ctStanModel or ctStanFit, and vector of values for free parameters.
#'
#' @param model either a ctStanFit or ctStanModel object.
#' @param parvalues vector of parameter values to assign to any free parameters in the model
#' @param timeinterval time interval to use for discrete time (dt) matrix calculations.
#'
#' @return A list containing various matrices related to a continuous time dynamic model. 
#' Matrices with "dt" in front refers to discrete time, "asym" refers to asymptotic (time interval = infinity), 
#' and "cor" refers to correlations. 
#' @export
#'
#' @examples
#' ctStanParMatrices(ctstantestfit,rep(0,16))
ctStanParMatrices <- function(model, parvalues, timeinterval=1){
  if(class(model)=='ctStanFit') model<- model$ctstanmodel
  if(class(model) !='ctStanModel') stop('not a ctStanModel')
  
  valuespec <- model$pars[is.na(model$pars$param),]
  parspec <- model$pars[!is.na(model$pars$param),]
  
  if(length(parvalues)!=nrow(parspec)) stop('length of parvalues != number of free params in model!')
  
  covchol2corchol<-function(mat, invert){ 
s=c()
    o=mat;
    
    for(i in 1:nrow(o)){ 
      for(j in min(i+1,nrow(mat)):nrow(mat)){
        o[j,i] = inv_logit(o[j,i])*2-1;  # can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=1; #change to adjust prior for correlations
    }
    if(invert==1) o = solve(o);
    
    for(i in 1:nrow(o)){
      s[i] = 1/sqrt(o[i,,drop=FALSE] %*% o[,i,drop=FALSE]);
      if(is.infinite(s[i])) s[i]=0;
    }
    o= diag(s) %*% o;
    return(o);
  }
  
sdcovchol2cov <- function(mat, cholesky){ 
    invert = 0; 
    if(nrow(mat) > 1){
      out=covchol2corchol(mat,invert); 
      out= diag(diag(mat)) %*% out
      }
    if(nrow(mat)==1) out = mat;
    
    if(cholesky==0) out = out %*% t(out);
    return(out);
  }
  
  
  
  for(mati in unique(model$pars$matrix)){
    assign(mati,matrix(NA, 
      nrow= max(model$pars$row[model$pars$matrix %in% mati]),
      ncol= max(model$pars$col[model$pars$matrix %in% mati]),
    ))
  }

  
  for(i in 1:nrow(valuespec)){
    eval(parse(text=paste0(valuespec$matrix[i], '[',valuespec$row[i],' ,', valuespec$col[i], '] <- ', valuespec$value[i])))
  }
  
  for(i in 1:nrow(parspec)){
    param <- parvalues[i]
    eval(parse(text=paste0(parspec$matrix[i], '[',parspec$row[i],' ,', parspec$col[i], '] <- ', parspec$transform[i])))
  }
  
DIFFUSION = sdcovchol2cov(DIFFUSION,0)
DIFFUSIONcor = suppressWarnings(stats::cov2cor(DIFFUSION))
DIFFUSIONcor[is.na(DIFFUSIONcor)] <- 0
T0VAR=sdcovchol2cov(T0VAR,0)
T0VARcor = suppressWarnings(stats::cov2cor(T0VAR))
T0VARcor[is.na(T0VARcor)] <- 0
MANIFESTVAR=MANIFESTVAR^2

DRIFTHATCH<-DRIFT %x% diag(nrow(DRIFT)) + diag(nrow(DRIFT)) %x% DRIFT
asymDIFFUSION<-matrix(-solve(DRIFTHATCH, c(DIFFUSION)), nrow=nrow(DRIFT))
asymDIFFUSIONcor = suppressWarnings(stats::cov2cor(asymDIFFUSION))
asymDIFFUSIONcor[is.na(asymDIFFUSIONcor)] <- 0

ln=model$latentNames
mn=model$manifestNames
tdn=model$TDpredNames

dimnames(DRIFT)=list(ln,ln)
dimnames(DIFFUSION)=list(ln,ln)
dimnames(asymDIFFUSION)=list(ln,ln)
rownames(CINT)=ln
rownames(MANIFESTMEANS)=mn
rownames(T0MEANS)=ln

dimnames(T0VAR)=list(ln,ln)
dimnames(asymDIFFUSION)=list(ln,ln)
dimnames(LAMBDA)=list(mn,ln)


dtDRIFT=expm(DRIFT * timeinterval)

dtDIFFUSION = asymDIFFUSION - (dtDRIFT %*% asymDIFFUSION %*% t(dtDRIFT ))
dtDIFFUSIONcor = cov2cor(dtDIFFUSION)

dtCINT = (solve(DRIFT) %*%(dtDRIFT - diag(nrow(DRIFT))) %*% (CINT))

asymCINT = -solve(DRIFT) %*% CINT

out<-list(DRIFT=DRIFT,dtDRIFT=dtDRIFT, T0VAR=T0VAR, T0VARcor=T0VARcor, 
  DIFFUSION=DIFFUSION, DIFFUSIONcor=DIFFUSIONcor, dtDIFFUSION=dtDIFFUSION, dtDIFFUSIONcor=dtDIFFUSIONcor,
  asymDIFFUSION=asymDIFFUSION, asymDIFFUSIONcor=asymDIFFUSIONcor, 
  CINT=CINT,dtCINT=dtCINT, asymCINT=asymCINT, T0MEANS=T0MEANS,
  MANIFESTMEANS=MANIFESTMEANS, LAMBDA=LAMBDA)

if('MANIFESTVAR' %in% model$pars$matrix) {
  dimnames(MANIFESTVAR)=list(mn,mn)
  out$MANIFESTVAR=MANIFESTVAR
  
}

if('TDPREDEFFECT' %in% model$pars$matrix) {
  dimnames(TDPREDEFFECT)=list(ln,tdn)
  out$TDPREDEFFECT<-TDPREDEFFECT
}
  
return(out)
}
