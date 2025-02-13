if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  
  if(F){  
    library(ctsem)
    library(testthat)
    library(ggplot2)
    library(data.table)
    library(gridExtra)
    set.seed(2)
    cores=5
    
    context("boothesscheck")
    
    test_that("boothesscheck1", {
      
      niter=500
      for(iteri in 1:niter){
        message(iteri)
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
        # tdat[tdat[,'id']==2,'TI1'] <- NA
        
        checkm<-suppressMessages(ctModel(type='ct',Tpoints=Tpoints,
          # MANIFESTVAR=diag(0.5,1),
          # DRIFT=matrix(c(-.3),nrow=1),
          # DIFFUSION=matrix(c(2),1),
          n.latent=n.latent,n.TDpred=n.TDpred,
          n.TIpred=n.TIpred,
          MANIFESTMEANS=matrix(0,nrow=n.manifest),
          CINT=matrix(c('cint1'),ncol=1),
          n.manifest=n.manifest,LAMBDA=diag(1)))
        
        trueval=c(as.numeric(gm$T0MEANS[1]),-.3,as.numeric(gm$DIFFUSION[1]),
          as.numeric(gm$MANIFESTVAR[1,1]),0,
          c(0,-1,0,0,5),
          c(1,0.6,0))
        
        # checkm$pars$indvarying <- FALSE
        
        # checkm$pars[c(-1,-7) ,c('TI1_effect')] <- FALSE
        
        tfit1<-ctStanFit(tdat,checkm,chains=1,cores=cores,optimize=TRUE,
          # plot=1,
          optimcontrol=list(bootstrapUncertainty=TRUE))
        s1=summary(tfit1,parmatrices=F)
        
        # test_isclose(s1$tipreds[2,'mean'],5,tol=.2)
        # test_isclose(s1$popsd[2,'mean'],.6,tol=.2)
        
        tfit2<-ctStanFit(tdat,checkm,optimize=TRUE,cores=cores,inits=tfit1$stanfit$rawest,
          optimcontrol=list(bootstrapUncertainty=F,stochastic=F))
        s2=summary(tfit2,parmatrices=F)
        
        # tfit3<-ctStanFit(tdat,checkm,chains=1,cores=cores,optimize=TRUE,
        #   optimcontrol=list(hessianType='stochastic',stochastic=F,stochasticHessianSamples=50,stochasticHessianEpsilon=1e-4))
        # s3=summary(tfit3,parmatrices=F)
        # s3$popmeans
        
        r1 <- data.frame(iteration=iteri,type='boot',trueval=trueval,
          rbind(s1$popmeans,s1$tipreds[,-6],s1$popsd,s1$rawpopcorr[,-6]))
        r2 <- data.frame(iteration=iteri,type='normal',trueval=trueval,
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
          # geom_point(aes(y=trueval),colour='black')+
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
      
      # tfit3<-ctStanFit(tdat,checkm,cores=cores,optimize=TRUE, inits=tfit1$stanfit$rawest,
      #   optimcontrol=list(hessianType='stochastic',stochasticHessianEpsilon=1e-1))
      # s3=summary(tfit3,parmatrices=F)
      # s3$popmeans
      
  })
}
}
