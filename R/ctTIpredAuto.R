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
  Tpoints=30
  n.manifest=1
  n.TDpred=0
  n.TIpred=1
  n.latent=1
  n.subjects=100
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
    ndat<-ctGenerate(gm,n.subjects=1,burnin=30,logdtsd=.4)
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
  
  fit1<-ctStanFit(tdat,checkm,chains=1,optimize=TRUE,cores=2,verbose=1,
    # intoverpop=F,plot=T,
    # savesubjectmatrices = F,plot=F,
    # init=init,
    optimcontrol=list(is=FALSE,stochastic=T,subsamplesize=1,carefulfit=F),
    nopriors=TRUE)
}

whichsubjectpars <- function(standata,subjects=NA){
  a1=standata$nparams+standata$nindvarying+
    (standata$nindvarying^2-standata$nindvarying)/2
  whichbase <- 1:a1
  if(standata$intoverpop ==0 && standata$nindvarying > 0){ #then there are subject pars
    whichsubjects <- a1+cseq(from=subjects,to=standata$nindvarying*standata$nsubjects,
      by=standata$nsubjects)
    whichbase <- c(whichbase,whichsubjects)
  }
  if(standata$ntipredeffects > 0) {
    tipredstart <- (a1+
        ifelse(standata$intoverpop,0,standata$nindvarying*standata$nsubjects)+1)
    whichbase <- c(whichbase,tipredstart:(tipredstart+standata$ntipredeffects+
        ifelse(standata$doonesubject >0,0,-1)))
  }
  return(whichbase)
}


scorecalc <- function(standata,est,stanmodel,subjectsonly=TRUE,
  returnsubjectlist=TRUE,cores=2){
  standata$dokalmanpriormodifier <- ifelse(subjectsonly, 1/standata$nsubjects,1/standata$ndatapoints)
  
  scores <- list()
  for(i in 1:standata$nsubjects){
    whichpars = whichsubjectpars(standata,i)
    scores[[i]]<-matrix(NA,length(whichpars),ifelse(subjectsonly,1,sum(standata$subject==i)))
    standata1 <- standatact_specificsubjects(standata,i)
    for(j in 1:ncol(scores[[i]])){
      standata1$llsinglerow=as.integer(ifelse(subjectsonly,0,j))
      sf <- stan_reinitsf(stanmodel,standata1,fast = TRUE)
      scores[[i]][,j] <- sf$grad_log_prob(
        upars=est[whichpars],
        adjust_transform = TRUE)
    }
  }
  
  if(subjectsonly) scores <- matrix(unlist(scores),nrow=length(scores[[i]]))
  if(!returnsubjectlist){ #return data.table
    if('list' %in% class(scores)){
      scores=lapply(scores,function(x) data.table(t(x)))
      scores=rbindlist(scores)
    } else scores <- t(scores)
   
  }
  return(scores)
}

ctTIauto <- function(fit,tipreds=NA){
  if(is.na(tipreds[1])) tipreds <- fit$standata$tipredsdata
  # colnames(tipreds) <- paste0('ti',1:ncol(tipreds))
  scores <- scorecalc(standata = fit$standata,
    est = fit$stanfit$rawest,stanmodel = fit$stanmodel)
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
    e=stan_constrainsamples(sm = fit$stanmodel,standata = fit$standata,
      samples = matrix(fit$stanfit$rawest,nrow=1),savescores=TRUE,quiet=TRUE,pcovn=2)
    p=sort(unique(fit$setup$matsetup$row[fit$setup$matsetup$indvarying>0]))# | fit$setup$matsetup$tipred]))
    firstsub <- rep(TRUE,fit$standata$ndatapoints) #which rows represent first rows per subject
    for(i in 2:fit$standata$ndatapoints){
      if(fit$standata$subject[i] == fit$standata$subject[i-1]) firstsub[i] <- FALSE
    }
    e$etasmooth <-  array(e$etaa[,3,,,drop=FALSE],dim=dim(e$etaa)[-2])
    states <- ctCollapse(e$etasmooth[,firstsub,p,drop=FALSE],1,mean)
    sc=list()
    for(i in 1:ncol(states)){
      sc[[i]]=summary(lm(states[,i] ~ tipreds))$coefficients[-1,,drop=FALSE]
    }
    for(i in 1:length(sc)){
      TIPREDEFFECTsetup[fit$setup$matsetup$param[fit$setup$matsetup$indvarying %in% i],] <- sc[[i]][,4]
    }
  }
  
  if(any(is.na(TIPREDEFFECTsetup))) warning('NA found, probably unused parameters?')
  TIPREDEFFECTsetup[is.na(TIPREDEFFECTsetup)] <- 1
  return(TIPREDEFFECTsetup)
}
ctIndVarAuto <- function(fit,aicthreshold = -2){
  if(requireNamespace('lme4',quietly=TRUE)){
    scores <- scorecalc(standata = fit$standata,est = fit$stanfit$rawest,
      stanmodel = fit$stanmodel,subjectsonly = FALSE)
    scores <- t(do.call(cbind,scores))
    scores <- scores[,1:fit$standata$nparams]
    colnames(scores) <- paste0('p',1:ncol(scores))
    sdat <- fit$standata
    sdat$savescores <- 1L
    # sdat$nopriors <- 1L
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
    try(suppressWarnings(suppressMessages({
      for(i in 1:ncol(scores)){
        # for(j in 1:ncol(tipreds)){
        # plot(sort(tipreds[,j]),scores[i,][order(tipreds[,j])],ylab=rownames(scores)[i],xlab=colnames(tipreds)[j])
        # }
        states <-c()
        dat=data.frame(scores=scores[,i],subject=fit$standata$subject,one=1,etaprior,etaprior2)
        f <- paste0('scores ~ (1|subject)')
        f1<-paste0('scores ~ (1|one)')
        f2<- paste0('+',c(colnames(etaprior),colnames(etaprior2)),collapse='+')
        l=lme4::lmer(data = dat,formula(paste0(f,f2)))
        s=summary(l)
        states<-rownames(s$coefficients[abs(s$coefficients[,3]) > 1.96 & 
            abs(s$coefficients[,1]) > .1,])[-1]
        f2 <- paste0(ifelse(length(states)>0,'+',''),states,collapse='+')
        l1 <- lme4::lmer(data = dat,formula(paste0(f,f2)))
        l2=lme4::lmer(data = dat,formula(paste0(f1,f2)),control=lme4::lmerControl(check.nlev.gtr.1="ignore"))
        a=as.matrix(anova(l,l1,l2))
        # a=a[2,]-a[1,]
        s1=summary(l1)
        if(s1$varcor$subject[1]^2 > .05 && a[2,'AIC']-a[1,'AIC'] < 3) indvarying[i] <- TRUE #if 5% or more variance explained 
        
        saic <-c()
        for(si in seq_along(states)){
          fs <- paste0(ifelse(length(states)>1,'+',''),states[-si],collapse='+')
          ls <- lme4::lmer(data = dat,formula(paste0(f,fs)))
          as <-  as.matrix(anova(l1,ls))
          saic[si] <- as[2,'AIC']-as[1,'AIC']
        }
        sumout <- s1$coefficients[,1,drop=FALSE]^2
        sumout <- cbind(sumout,c(0,saic))
        sumout <- rbind(sumout, c(s1$varcor$subject[1]^2,a[1,'AIC']-a[2,'AIC']))[-1,,drop=FALSE]
        rownames(sumout)[nrow(sumout)] <- 'random'
        colnames(sumout)[2] <- 'AICdiff'
        out[[i]] <- sumout
        statelist[[i]] <- states
      }
    })),silent=TRUE)
    
    names(out)[1:fit$standata$nparams]<-fit$setup$matsetup$parname[match(1:fit$standata$nparams,fit$setup$matsetup$param)]
    o=sapply(out,min)
    out <- out[order(o)]
    out <- lapply(out,function(x) x[x[,2]< aicthreshold,,drop=FALSE])
    out <- out[lapply(out,length)>0]
    return(out)
  } else stop('lme4 package needed -- install.packages("lme4")')
}
