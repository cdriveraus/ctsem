checkTIauto <- function(){
  Tpoints=20
  n.manifest=1
  n.TDpred=0
  n.TIpred=1
  n.latent=1
  n.subjects=80
  TI1 <- rnorm(n.subjects)
  gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
    n.TDpred=n.TDpred,n.manifest=n.manifest,
    MANIFESTVAR=diag(0.5,1),
    LAMBDA=diag(1,1),
    DRIFT=matrix(c(-.3),nrow=1),
    DIFFUSION=matrix(c(2),1),
    T0VAR=diag(10,1))
  
  for(i in 1:n.subjects){
    gm$CINT[1,1] <- TI1[i]*5+rnorm(1,0,.6)
    ndat<-ctGenerate(gm,n.subjects=1,burnin=30,wide=FALSE,logdtsd=.4)
    ndat <- cbind(ndat,TI1[i])
    ndat[,1] <- i
    if(i>1) tdat <- rbind(tdat,ndat) else tdat <- ndat
  }
  colnames(tdat)[4] <- 'TI1'
  
  tdat[2,'Y1'] <- NA
  tdat[tdat[,'id']==2,'TI1'] <- NA
  
  checkm<-ctModel(type='stanct',Tpoints=Tpoints,
    MANIFESTVAR=diag(0.5,1),
    DRIFT=matrix(c(-.3),nrow=1),
    DIFFUSION=matrix(c(2),1),
    n.latent=n.latent,n.TDpred=n.TDpred,
    n.TIpred=n.TIpred,
    MANIFESTMEANS=matrix(0,nrow=n.manifest),
    CINT=matrix(c('cint1'),ncol=1),
    n.manifest=n.manifest,LAMBDA=diag(1))
  
  # checkm$pars$indvarying <- FALSE
  
  checkm$TIpredAuto <- 1L
  
  tfit2<-ctStanFit(tdat,checkm,chains=1,optimize=TRUE,cores=1,verbose=1,
    optimcontrol=list(is=FALSE),nopriors=TRUE,
    nlcontrol=list(nldynamics=TRUE,nlmeasurement=TRUE))
}


scorecalc <- function(fit,subjectsonly=TRUE){
  fit$standata$nopriors=1L
  
  scores <- list()
  for(i in 1:fit$standata$nsubjects){
    scores[[i]]<-matrix(NA,length(fit$stanfit$rawest),ifelse(subjectsonly,1,sum(fit$standata$subject==i)))
    standata <- standatact_specificsubjects(fit$standata,i)
    for(j in 1:ncol(scores[[i]])){
      standata$llsinglerow=as.integer(ifelse(subjectsonly,0,j))
    sf <- stan_reinitsf(fit$stanmodel,standata,fast = TRUE)
    scores[[i]][,j] <- sf$grad_log_prob(upars=fit$stanfit$rawest,adjust_transform = TRUE)
    }
  }
  # browser()
  if(subjectsonly) scores <- matrix(unlist(scores),nrow=length(scores[[i]]))
  return(scores)
}

ctTIauto <- function(fit,tipreds=NA){
  if(is.na(tipreds[1])) tipreds <- fit$standata$tipredsdata
  # colnames(tipreds) <- paste0('ti',1:ncol(tipreds))
  scores <- scorecalc(fit)
  scores <- scores[1:fit$standata$nparams,,drop=FALSE]
  rownames(scores) <- paste0('p',1:nrow(scores))
  # matchindex <- match(1:fit$standata$nparams,fit$setup$matsetup$param)
  # rownames(scores)[1:fit$standata$nparams] <- fit$standata$matsetup$parname[match(1:fit$standata$nparams,fit$setup$matsetup$param)]
  sc <- list()
  for(i in 1:nrow(scores)){
    # for(j in 1:ncol(tipreds)){
    # plot(sort(tipreds[,j]),scores[i,][order(tipreds[,j])],ylab=rownames(scores)[i],xlab=colnames(tipreds)[j])
    # }
    sc[[i]]=summary(lm(scores[i,] ~ tipreds))$coefficients
  }
  names(sc)[1:fit$standata$nparams]<-paste0('p',1:nrow(scores)) #fit$setup$matsetup$parname[match(1:fit$standata$nparams,fit$setup$matsetup$param)]
  
  s2=lapply(sc,function(x) {
    x=x[-1,,drop=FALSE]
    rownames(x) <- gsub('^tipreds','',rownames(x))
    rownames(x) <- paste0('ti',1:nrow(x))
    return(x)
  })
  
  TIPREDEFFECTsetup = matrix(NA,length(s2),nrow(s2[[1]]))
  for(i in 1:length(s2)){
    TIPREDEFFECTsetup[i,] <- s2[[i]][,4]
  }
  
  if(fit$standata$nindvarying > 0 && fit$standata$intoverpop > 0){
    fit$stanfit$rawposterior <- matrix(fit$stanfit$rawest,1)
    fit$setup$matsetup <- data.frame(fit$standata$matsetup)
    fit$standata$savescores <- 1L
    e=stan_constrainsamples(sm = fit$stanmodel,standata = fit$standata,samples = matrix(fit$stanfit$rawest,nrow=1))
    p=sort(unique(fit$setup$matsetup$row[fit$setup$matsetup$indvarying>0]))# | fit$setup$matsetup$tipred]))
    firstsub <- rep(TRUE,fit$standata$ndatapoints) #which rows represent first rows per subject
    for(i in 2:fit$standata$ndatapoints){
      if(fit$standata$subject[i] == fit$standata$subject[i-1]) firstsub[i] <- FALSE
    }
    states <- ctCollapse(e$etasmooth[,firstsub,p,drop=FALSE],1,mean)
    sc=list()
    for(i in 1:ncol(states)){
      sc[[i]]=summary(lm(states[,i] ~ tipreds))$coefficients[-1,,drop=FALSE]
    }
    for(i in 1:length(sc)){
      TIPREDEFFECTsetup[fit$setup$matsetup$param[fit$setup$matsetup$indvarying %in% i],] <- sc[[i]][,4]
    }
  }
  
  # s2=lapply(s2,function(x) x=x[x[,4] < .05,,drop=FALSE])
  # s2=s2[which(unlist(lapply(s2,function(x) length(x) > 0)))]
  return(TIPREDEFFECTsetup)
}

# parallelStanSetup <- function(cl, standata,split=TRUE){
#   cores <- length(cl)
#   stansubjectindices <- split(unique(standata$subject),sort(unique(standata$subject) %% min(standata$nsubjects,cores)))
#   
#   parallel::clusterExport(cl,c('standata'),envir = environment())
#   parallel::clusterApply(cl,stansubjectindices,function(subindices) {
#     # require(Rcpp)
#     library(ctsem)
#     if(length(subindices) < length(unique(standata$subject))) standata <- standatact_specificsubjects(standata,subindices)
#     if(!1 %in% subindices) standata$nopriors <- 1L
#     g = eval(parse(text=paste0('gl','obalenv()'))) #avoid spurious cran check -- assigning to global environment only on created parallel workers.
#     assign('smf',stan_reinitsf(sm,standata),pos = g)
#     NULL
#   })
#   NULL
# }
