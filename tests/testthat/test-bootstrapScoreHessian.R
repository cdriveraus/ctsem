if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  
  if(F){  

    
    context("boothesscheck")
    
    test_that("boothesscheck1", {
      library(ctsem)
      library(testthat)
      library(ggplot2)
      library(data.table)
      library(gridExtra)
      set.seed(2)
      cores=10
      niter=500
      
      #generating model / simulated data / true pars
      Tpoints=30
      n.manifest=1
      n.TDpred=0
      n.TIpred=1
      n.latent=1
      n.subjects=50
      TI1 <- rbinom(n.subjects,1,.5)
      gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
        n.TDpred=n.TDpred,n.manifest=n.manifest,
        MANIFESTVAR=diag(0.5,1),
        LAMBDA=diag(1,1),T0MEANS=10,
        DRIFT=matrix(c(-.3),nrow=1),
        DIFFUSION=matrix(c(2),1),
        T0VAR=diag(1,1))
      
      for(i in 1:n.subjects){
        gm$CINT[1,1] <- TI1[i]*5+rnorm(1,0,.6)
        gm$DRIFT[1,1] <- -.3+TI1[i]*-1
        ndat<-suppressMessages(ctGenerate(gm,n.subjects=1,burnin=0,logdtsd=.4,dtmean = .2))
        ndat <- cbind(ndat,TI1[i])
        ndat[,1] <- i
        if(i>1) tdat <- rbind(tdat,ndat) else tdat <- ndat
      }
      colnames(tdat)[4] <- 'TI1'
      
      tdat[2,'Y1'] <- NA
      
      checkm<-suppressMessages(ctModel(type='ct',Tpoints=Tpoints,
        n.latent=n.latent,n.TDpred=n.TDpred,
        n.TIpred=n.TIpred,
        MANIFESTMEANS=matrix(0,nrow=n.manifest),
        CINT=matrix(c('cint1'),ncol=1),
        n.manifest=n.manifest,LAMBDA=diag(1)))

      
      tfit1<-ctStanFit(tdat,checkm,cores=cores,optimize=TRUE,
        optimcontrol=list(bootstrapUncertainty=TRUE))
      
      #use the estimates as true pars going forward
      truepars <- tfit1$stanfit$rawest
      
      tfit1 <- ctStanGenerateFromFit(fit = tfit1, nsamples = niter, cores=cores)
      
      
      for(iteri in 1:niter){
        message(iteri)
       
        y <- array(tfit1$generated$Y[iteri,,],dim=dim(tfit1$generated$Y[1,,,drop=FALSE])[-1])
        colnames(y) <- checkm$manifestNames
        dat <- tdat
        dat[,checkm$manifestNames] <- y
        
        fitboot<-ctStanFit(tdat,checkm,cores=cores,
          optimcontrol=list(bootstrapUncertainty=TRUE))
        sboot=summary(fitboot,parmatrices=F)

        
        fithess<-ctStanFit(tdat,checkm,optimize=TRUE,cores=cores,inits=tfit1$stanfit$rawest,
          optimcontrol=list(bootstrapUncertainty=F,stochastic=F))
        s2=summary(tfit2,parmatrices=F)

        r1 <- data.frame(iteration=iteri,type='boot',truepars=truepars,
          rbind(s1$popmeans,s1$tipreds[,-6],s1$popsd,s1$rawpopcorr[,-6]))
        r2 <- data.frame(iteration=iteri,type='normal',truepars=truepars,
          rbind(s2$popmeans,s2$tipreds[,-6],s2$popsd,s2$rawpopcorr[,-6]))
        
        r1$par <- rownames(r1)
        r2$par <- rownames(r2)
        
        resiter <- rbind(r1,r2)
        
        if(iteri==1) res <- resiter else res <- rbind(res,resiter)
        
        res2=data.table(res)
        res2[,avgMean:=mean(mean),by='par']
        res2[,coverage:=mean(X2.5.<=trueval & X97.5.>=trueval),by=c('type','par')]
        
        g1=ggplot(res2,aes(x=par,y=mean-trueval,colour=type))+
            geom_boxplot(outliers = F)+ 
          geom_boxplot(aes(y=X2.5.-trueval),outliers = F,linetype='dotted')+
          geom_boxplot(aes(y=X97.5.-trueval),outliers = F,linetype='dotted')+
            theme_bw()+
            geom_hline(yintercept=0,linetype='dashed')+
            theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position='bottom')
        
        g2=ggplot(res2[iteration==1,],aes(x=par,y=coverage,colour=type))+
          geom_point(alpha=.7)+ #rotate x axis text
            theme_bw()+
          geom_hline(yintercept=.95,linetype='dashed')+
          theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position='bottom')
        
        #plot g1 next to g2
        print(grid.arrange(g1,g2,nrow=2))
          
      }

      
  })
}
}
