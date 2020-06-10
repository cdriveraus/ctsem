if(1==99){
set.seed(1)
    Tpoints=50
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
      gm$CINT[1,1] <- TI1[i]*.5+rnorm(1,.4,1.6)
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
      DRIFT=matrix(c(-.3),nrow=1),
      DIFFUSION=matrix(c(2),1),
      n.latent=n.latent,n.TDpred=n.TDpred,
      n.TIpred=n.TIpred,
      MANIFESTMEANS=matrix(0,nrow=n.manifest),
      CINT=matrix(c('cint1'),ncol=1),
      n.manifest=n.manifest,LAMBDA=diag(1))
    
    # checkm$pars$indvarying <- FALSE
    
    checkm$pars[c(-1,-7) ,c('TI1_effect')] <- FALSE
    

    
    tfit2<-ctStanFit(tdat,checkm,chains=1,optimize=TRUE,cores=6,verbose=1,
      optimcontrol=list(is=FALSE),nopriors=TRUE,
      nlcontrol=list(nldynamics=TRUE,nlmeasurement=TRUE))
    
    summary(tfit2)
    
    sm1<-tfit2$stanmodel
    sd1=tfit2$standata
    sd1$verbose=2L
    sf1<-stan_reinitsf(sm1,sd1)
    e1 <- ctExtract(tfit2)
    p1=tfit2$stanfit$rawest
    
    sm2<-tfit2$stanmodel
    sd2=tfit2$standata
    sd2$verbose=2L
    sf2<-stan_reinitsf(sm2,sd2)
    e2 <- ctExtract(tfit2)
    p2=tfit2$stanfit$rawest
    
    sink('sink1.txt')
    rstan::log_prob(sf1,p1)
    sink()
    sink('sink4.txt')
    rstan::log_prob(sf2,p2)
    sink()
    
    p1=c(0.13676877,  0.03638180 , 0.39071061 ,-0.37734559 , 1.51513814,  0.13338917,  0.04727474)
    
    sc<-list()
    for(i in 1:nrow(scores)){
      # for(j in 1:ncol(tipreds)){
      # plot(sort(tipreds[,j]),scores[i,][order(tipreds[,j])],ylab=rownames(scores)[i],xlab=colnames(tipreds)[j])
      # }
      dat=data.frame(scores=scores[i,],subject=fit$standata$subject,one=1)
      l=lmer(data = dat,scores ~ (1|subject))
      l2=lmer(data = dat,scores ~ (1|one),control=lmerControl(check.nlev.gtr.1="ignore"))
      message(fit$setup$matsetup$parname[fit$setup$matsetup$param %in% i],' = ',diff(anova(l2,l)$AIC))
    }
}
