#' Returns model equation and other matrices from a ctStanFit object, and vector of values for free parameters.
#'
#' @param fit ctStanFit object.
#' @param parvalues vector of parameter values to assign to free parameters in the model
#' @param timeinterval time interval to use for discrete time (dt) matrix calculations.
#' @param sf stanfit object. Generally not necessary, but for repeated calls to this function, can speed things up.
#'
#' @return A list containing various matrices related to a continuous time dynamic model. 
#' Matrices with "dt" in front refers to discrete time, "asym" refers to asymptotic (time interval = infinity), 
#' and "cor" refers to correlations. 
#' @export
#'
#' @examples
#' ctStanParMatrices(ctstantestfit,rep(0,16))
ctStanParMatrices <- function(fit, parvalues, timeinterval=1, sf=NA){

  if(class(fit) !='ctStanFit') stop('not a ctStanFit object')
  model <- fit$ctstanmodel

  if(length(parvalues)!=fit$data$nparams) stop('length of parvalues != number of free params (',fit$data$nparams,') in model!')
  
#   covchol2corchol<-function(mat, invert){ 
# s=c()
#     o=mat;
#     
#     for(i in 1:nrow(o)){ 
#       for(j in min(i+1,nrow(mat)):nrow(mat)){
#         o[j,i] = inv_logit(o[j,i])*2-1;  # can change cor prior here
#         o[i,j] = o[j,i];
#       }
#       o[i,i]=1; #change to adjust prior for correlations
#     }
#     if(invert==1) o = solve(o);
#     
#     for(i in 1:nrow(o)){
#       s[i] = 1/sqrt(o[i,,drop=FALSE] %*% o[,i,drop=FALSE]);
#       if(is.infinite(s[i])) s[i]=0;
#     }
#     o= diag(s) %*% o;
#     return(o);
#   }
#   
# sdcovchol2cov <- function(mat, cholesky){ 
#     invert = 0; 
#     if(nrow(mat) > 1){
#       out=covchol2corchol(mat,invert); 
#       out= diag(diag(mat)) %*% out
#       }
#     if(nrow(mat)==1) out = mat;
#     
#     if(cholesky==0) out = out %*% t(out);
#     return(out);
#   }
  
  if(all(is.na(sf[1]))) sf <- fit$stanfit
  
npars <- try(get_num_upars(sf),silent = TRUE) #$stanmodel)

if(class(npars)=='try-error'){ #in case R has been restarted or similar
  standataout <- fit$data
  standataout<-unlist(standataout)
  standataout[is.na(standataout)] <- 99999
  standataout <- utils::relist(standataout,skeleton=fit$data)
  suppressOutput(sf <- suppressWarnings(sampling(fit$stanmodel,data=standataout,iter=1,control=list(max_treedepth=1),chains=1)))
  npars <- get_num_upars(sf)
  pars <- c(parvalues,rep(0,npars - fit$data$nparams))
  sf <- constrain_pars(sf, pars)
} else { #if no problem getting npars...
  pars <- c(parvalues,rep(0,npars - fit$data$nparams))
  sf <- try(constrain_pars(sf, pars))
}

for(m in fit$setup$basematrices){
  assign(m,ctCollapse(sf[[paste0('pop_',m)]],1,mean))
}

choltrue <- !as.logical(fit$data$lineardynamics * fit$data$intoverstates)

if(choltrue) DIFFUSION = msquare(DIFFUSION) #sdcovchol2cov(DIFFUSION,0)
DIFFUSIONcor = suppressWarnings(stats::cov2cor(DIFFUSION))
DIFFUSIONcor[is.na(DIFFUSIONcor)] <- 0



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



dimnames(asymDIFFUSION)=list(ln,ln)
dimnames(LAMBDA)=list(mn,ln)


dtDRIFT=expm(DRIFT * timeinterval)

dtDIFFUSION = asymDIFFUSION - (dtDRIFT %*% asymDIFFUSION %*% t(dtDRIFT ))
dtDIFFUSIONcor = cov2cor(dtDIFFUSION)

dtCINT = (solve(DRIFT) %*%(dtDRIFT - diag(nrow(DRIFT))) %*% (CINT))

asymCINT = -solve(DRIFT) %*% CINT



if(choltrue) T0VAR=msquare(T0VAR)
T0VARcor = suppressWarnings(stats::cov2cor(T0VAR))
T0VARcor[is.na(T0VARcor)] <- 0
dimnames(T0VAR)=list(ln,ln)

rownames(T0MEANS)=ln

# for(i in row(statspec)){
#   if(statspec$matrix[i] =='T0VAR') {
#     eval(parse(text=paste0(statspec$matrix[i], '[',statspec$row[i],' ,', statspec$col[i], '] <- ', 
#     'asymDIFFUSION[',statspec$row[i],' ,', statspec$col[i], ']')))
#     eval(parse(text=paste0(statspec$matrix[i], 'cor[',statspec$row[i],' ,', statspec$col[i], '] <- ', 
#       'asymDIFFUSIONcor[',statspec$row[i],' ,', statspec$col[i], ']')))
#   }
#   if(statspec$matrix[i] =='T0MEANS') eval(parse(text=paste0(statspec$matrix[i], '[',statspec$row[i],' ,', statspec$col[i], '] <- ', 
#     'asymCINT[',statspec$row[i],' ,', statspec$col[i], ']')))
# }

T0VAR[upper.tri(T0VAR)] = t(T0VAR)[upper.tri(T0VAR)]
T0VARcor[upper.tri(T0VAR)] = t(T0VARcor)[upper.tri(T0VAR)]


out<-list(DRIFT=DRIFT,dtDRIFT=dtDRIFT, T0VAR=T0VAR, T0VARcor=T0VARcor, 
  DIFFUSION=DIFFUSION, DIFFUSIONcor=DIFFUSIONcor, dtDIFFUSION=dtDIFFUSION, dtDIFFUSIONcor=dtDIFFUSIONcor,
  asymDIFFUSION=asymDIFFUSION, asymDIFFUSIONcor=asymDIFFUSIONcor, 
  CINT=CINT,dtCINT=dtCINT, asymCINT=asymCINT, T0MEANS=T0MEANS,
  MANIFESTMEANS=MANIFESTMEANS, LAMBDA=LAMBDA)


 if(choltrue) MANIFESTVAR=msquare(MANIFESTVAR)
  dimnames(MANIFESTVAR)=list(mn,mn)
  out$MANIFESTVAR=MANIFESTVAR


if('TDPREDEFFECT' %in% model$pars$matrix) {
  dimnames(TDPREDEFFECT)=list(ln,tdn)
  out$TDPREDEFFECT<-TDPREDEFFECT
}
  
return(out)
}
