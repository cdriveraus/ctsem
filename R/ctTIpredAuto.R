ctStanParnamesRaw <- function(fit){
  if(!class(fit) %in% 'ctStanFit' && !class(fit$stanfit) %in% 'list') stop('Not an optimized ctStanFit model!')
  ms <- fit$setup$matsetup
  ms <- ms[ms$param > 0 & ms$when %in% 0 & ms$copyrow < 1,]
  ms <- ms[order(ms$param),]
  rawpopmeans<-ms$parname
  iv <- ms$parname[ms$indvarying >0]
  rawpopsd <-c()
  rawpopucorr <- c()
  if(length(iv) > 0){
    rawpopsd <- paste0('rawpopsd_',iv)
    rawpopucorr <- matrix(paste0('rawpopucorr_',rep(iv,times=length(iv)),'_',rep(iv,each=length(iv))),length(iv),length(iv))
    rawpopucorr <- rawpopucorr[lower.tri(rawpopucorr)]
  } 
  tipredeffects <- c()
  if(fit$standata$ntipred > 0){
    tipredeffects <- rep(NA,max(fit$standata$TIPREDEFFECTsetup))
    counter <- 0
    for(j in 1:ncol(fit$standata$TIPREDEFFECTsetup)){
      for(i in 1:nrow(fit$standata$TIPREDEFFECTsetup)){
        if(fit$standata$TIPREDEFFECTsetup[i,j] > 0){
          counter <- counter + 1
          tipredeffects[counter] <- paste0(rawpopmeans[i],'_', fit$ctstanmodel$TIpredNames[j])
        }
      }
    }
  }
  return(list(rawpopmeans=rawpopmeans,rawpopsd=rawpopsd,rawpopucorr=rawpopucorr,tipredeffects=tipredeffects))
}

TIPREDEFFECTnames <- function(fit){
  if(fit$standata$ntipred > 0){
    ms <- fit$setup$matsetup
    tie <- fit$standata$TIPREDEFFECTsetup
    colnames(tie) <- paste0(fit$ctstanmodel$TIpredNames)
    pars <- ms[ms$param > 0 & ms$when %in% 0 & ms$copyrow < 1,]
    pars <- pars[order(pars$param),]
    rownames(tie) <- pars$parname
    return(tie)
  } else return(c())
}

popcovnames <- function(fit){
  if(fit$standata$nindvarying > 0){
    ms <- fit$setup$matsetup
    pars <- ms[ms$param > 0 & ms$when %in% 0 & ms$copyrow < 1 & ms$indvarying > 0,]
    pars <- pars[order(pars$param),]
    popcov <- matrix(NA,fit$standata$nindvarying ,fit$standata$nindvarying,
      dimnames = list( pars$parname, pars$parname) )
    return(popcov)
  } else return(c())
}




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
    gm$CINT[1,1] <- TI1[i]*.5+rnorm(1,0,.6)
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
    # DRIFT=matrix(c(-.3),nrow=1),
    # DIFFUSION=matrix(c(2),1),
    n.latent=n.latent,n.TDpred=n.TDpred,
    n.TIpred=n.TIpred,
    MANIFESTMEANS=matrix(0,nrow=n.manifest),
    CINT=matrix(c('cint1'),ncol=1),
    n.manifest=n.manifest,LAMBDA=diag(1))
  
  checkm$pars$indvarying[!checkm$pars$matrix %in% 'T0MEANS'] <- FALSE
  
  checkm$TIpredAuto <- 1L
  
  fit<-ctStanFit(tdat,checkm,chains=1,optimize=TRUE,cores=1,verbose=1,
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
ctIndVarAuto <- function(fit){
  scores <- scorecalc(fit,subjectsonly = FALSE)
  scores <- t(do.call(cbind,scores))
  scores <- scores[,1:fit$standata$nparams]
  colnames(scores) <- paste0('p',1:ncol(scores))
  sdat <- fit$standata
  sdat$savescores <- 1L
  sdat$nopriors <- 1L
  sdat$popcovn <- 5L
  e=rstan::constrain_pars(stan_reinitsf(fit$stanmodel,sdat),fit$stanfit$rawest)
  colnames(e$etaprior) <- paste0('etaprior',1:ncol(e$etaprior))
  etaprior <- scale(e$etaprior[,1:fit$ctstanmodel$n.latent,drop=FALSE])
  etaprior2 <- scale(etaprior^2)
  colnames(etaprior2) <- paste0(colnames(etaprior),'_sq')
  scores <- scale(scores)
  statelist <- list()
  indvarying <- c()
  out <-list()
  for(i in 1:ncol(scores)){
    # for(j in 1:ncol(tipreds)){
    # plot(sort(tipreds[,j]),scores[i,][order(tipreds[,j])],ylab=rownames(scores)[i],xlab=colnames(tipreds)[j])
    # }
    states <-c()
    dat=data.frame(scores=scores[,i],subject=fit$standata$subject,one=1,etaprior,etaprior2)
    f <- paste0('scores ~ (1|subject)')
    f1<-paste0('scores ~ (1|one)')
    f2<- paste0('+',c(colnames(etaprior),colnames(etaprior2)),collapse='+')
    l=lmer(data = dat,formula(paste0(f,f2)))
    s=summary(l)
    states<-rownames(s$coefficients[abs(s$coefficients[,3]) > 1.96 & 
        abs(s$coefficients[,1]) > .1,])[-1]
    f2 <- paste0(ifelse(length(states)>0,'+',''),states,collapse='+')
    l1 <- lmer(data = dat,formula(paste0(f,f2)))
    l2=lmer(data = dat,formula(paste0(f1,f2)),control=lmerControl(check.nlev.gtr.1="ignore"))
    a=as.matrix(anova(l,l1,l2))
    # a=a[2,]-a[1,]
    s1=summary(l1)
    if(s1$varcor$subject[1]^2 > .05 && a[2,'AIC']-a[1,'AIC'] < 3) indvarying[i] <- TRUE #if 5% or more variance explained 
    
    saic <-c()
    for(si in seq_along(states)){
      fs <- paste0(ifelse(length(states)>1,'+',''),states[-si],collapse='+')
      ls <- lmer(data = dat,formula(paste0(f,fs)))
      as <-  as.matrix(anova(l1,ls))
      saic[si] <- as[1,'AIC']-as[2,'AIC']
    }
    sumout <- s1$coefficients[,1,drop=FALSE]^2
    sumout <- cbind(sumout,c(0,saic))
    sumout <- rbind(sumout, c(s1$varcor$subject[1]^2,a[1,'AIC']-a[2,'AIC']))[-1,,drop=FALSE]
    rownames(sumout)[nrow(sumout)] <- 'random'
    colnames(sumout)[2] <- 'AICdiff'
    out[[i]] <- sumout
    statelist[[i]] <- states
  }
  names(out)[1:fit$standata$nparams]<-fit$setup$matsetup$parname[match(1:fit$standata$nparams,fit$setup$matsetup$param)]
  return(out)
}
