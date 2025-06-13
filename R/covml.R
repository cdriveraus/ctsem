
# if(1==99){ #random cor tests
# set.seed(1)
#   d=6
#   rm=matrix(rnorm(d^2,0,1),d)
#   # rm[diag(d)==1] <- 10
#   rm = rm%*% t(rm)
#   cholrm=t(chol(rm))
#   ndat=(t(cholrm %*% matrix(rnorm(1000*d,0,1),d)))
#   ndat[sample(1:length(ndat),ceiling(.1*length(ndat)),replace=FALSE)] <- NA
#   # ndat[-1:-2,c(3,4)] <- NA
#   ecov=cov(ndat,use="pairwise.complete.obs")
# 
#   system.time({ctfit=ctsem:::covml(ndat,reg = 0,independent = F,
#     corpriortype = 0,hmc=F,verbose = 0)})
#   ctcov=ctfit$estimate$covm
#   round(ctcov - ecov,2)
# 
#   system.time({mxfit=OpenMx::mxRefModels(data.frame(ndat),run=TRUE)})
#   mxcov <- tcrossprod(mxfit$Saturated$ltCov$values)
#   round(mxcov - ecov,3)
# 
#   #compare ctsem to openmx saturated
#   round(ctcov - mxcov,2)
#   round(cov2cor(ctcov) - cov2cor(mxcov),2)
# 
#   #compare ctsem to openmx variances
#   plot(1:d,
#     diag(ctcov) - diag(ecov),
#     type='b')
#   points(1:d,diag(mxcov) - diag(ecov),type='b',col='red')
# 
# rm
# ctfit$lp - mxfit$Saturated$output$Minus2LogLikelihood*-.5
## need to get working with missings:
# sum(mvtnorm::dmvnorm(x = ndat,mean = ctfit$estimate$mu,sigma = ctcov,log = TRUE))

## old function
# constraincorsqrt=function(rawcor, d){
#   o=matrix(NA,d,d)
#   counter=0;
#   for(i in 1:d){
#     for(j in 1:d){
#       if(j > i){
#         counter=counter+1;
#         o[j,i] =  inv_logit(rawcor[counter])*2-1;  #// can change cor prior here
#         o[i,j] = o[j,i];
#       }
#     }
#     o[i,i]=1;
#     o[i,] = o[i,] / sqrt(sum((o[i,])^2)+1e-10);
#   }
#   return(o);
# }
# 
# sdvec=exp(rnorm(6,0,1))
# corsqrt = constraincorsqrt(rnorm((6^2-6)/2,0,1),6)
# corm = tcrossprod(corsqrt)
# covm=diag(sdvec) %*% corm %*% diag(sdvec)
# round(cov2cor(covm)-corm,3)
# covsqrt = diag(sdvec) %*% corsqrt
# round(cov2cor(tcrossprod(covsqrt ))-corm)
# 
# }

covdata <- function(ndat,reg,independent=FALSE,corpriortype=1L){
  sdat = as.matrix(ndat)
  d=ncol(sdat)
  n=nrow(sdat)
  obs = matrix(!is.na(sdat),ncol=d)
  nobs=array(as.integer(apply(obs,1,sum)),dim=c(n))
  sdat[!obs] <- 0
  obs=t(apply(obs,1, function(x) c(which(x),rep(0,d-sum(x)))))
  
  covdata <- list(d=as.integer(d),n=as.integer(n),
    dat=sdat,
    corpriortype=as.integer(corpriortype),
    reg=(reg),
    indep=as.integer(independent),
    obs=array(as.integer(obs),dim=c(n,d)),
    nobs=nobs,
    symm=0L)
}

covml <- function(dat,reg=0,verbose=0,hmc=FALSE,
  independent=FALSE,corpriortype=2L,tol=1e-8){
  
  covdata=covdata(dat,reg,independent,corpriortype)
  d=covdata$d
  scovf <- suppressMessages(sampling(object = stanmodels$cov,iter=1,chains=0,check_data=FALSE,
    data=covdata))
  
  target <- function(pars){
    lp=try(rstan::log_prob(scovf,pars, gradient=TRUE))
    if(is.nan(lp) || 'try-error' %in% class(lp) ) {
      lp <- -999999999
      attributes(lp)$gradient <- rnorm(length(pars))
    }
    if(any(is.na(attributes(lp)$gradient))){
      message('NA gradient encountered...')
      attributes(lp)$gradient[is.na(attributes(lp)$gradient)] <- 0
    }
    if(verbose > 0 || verbose==TRUE) print(lp[1],digits = 20)
    # plot(pars)
    return(lp)
  }
  
  mizelist <- list(fn= function(x) -target(x)[1], 
    gr=function(x) -attributes(target(x))$gradient,
    fg=function(x){
      o=target(x)
  out=list(fn=-o[1],gr=-attributes(o)$gradient)
  return(out)
    }
  )
  
  init=c(apply(as.matrix(dat),2,mean,na.rm=TRUE),
    log(c(apply(as.matrix(dat),2, sd,na.rm=TRUE))))
  if(!independent) init=c(init,rnorm((d^2-d)/2,0,.01))
  init[is.na(init) | is.infinite(init)]=0

  if(!hmc){
    # covfit=sgd(init = init,fitfunc = target,plot = 1)
  covfit=mize(par = init,fg=mizelist,memory=20,max_iter=10000,
  #   # line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
    abs_tol=tol,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
  } else{
    # browser()
    covfit=stanWplot(object = stanmodels$cov,iter=2000,chains=4,cores=4,check_data=FALSE,
      data=covdata,init_r=.01)
    # covfit=sampling(object = stanmodels$cov,iter=2000,chains=2,cores=2,check_data=FALSE,
    #   data=covdata)
    
  }

  cp=rstan::constrain_pars(object = scovf,covfit$par)
  dimnames(cp$covm)=list(colnames(dat),colnames(dat))
  
  return(list(dat=covdata,fit=covfit, estimate=cp,
    ll=sum(cp$llrow,na.rm=TRUE), 
    lp=rstan::log_prob(scovf,covfit$par)))
}
