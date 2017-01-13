# library(ctsem)
# library(testthat)
# library(rstan)
# 
# context("Bayes vs Freq")
# 
# 
# test_that("time calc", {
#   
  
  pars=cbind(c('drift11',-.1),
    c('tdpredeffect11', .1),
    c('tdpredeffect12', -3),
    c('diffusion11',.3),
    c('manifestvar11',.2),
    c('t0var11',3),
    c('cint1',4),
    c('manifestmeans2',3),
    c('lambda32',-2)
  )
  out=c()
  
  for(pari in 1:ncol(pars)){
    for(parvali in 2:sum(!is.na(pars[,pari]))){
      for(nsubjects in c(5,20)){
        for(tpoints in c(6,25)){
          
          for(i in 1:ncol(pars)){
            assign(pars[1,i],pars[2,i])
          }
          
          
          
          varpar=pari %in% c(4,5,6)
          
          nlatent=2
          nmanifest=3
          ntdpred=3
          TDPREDMEANS=matrix(0,ntdpred*(tpoints),1)
          TDPREDMEANS[runif(10,0,tpoints*ntdpred)]=1

          
          
          genm=ctModel(Tpoints=tpoints,
            n.latent=nlatent, n.manifest=nmanifest, n.TDpred=ntdpred,
            LAMBDA=matrix(c(1, 0,0,1,.6,lambda32), nrow=nmanifest, ncol=nlatent),
            DRIFT=matrix(c(drift11, 1, 0, -0.5), byrow=TRUE, nrow=nlatent, ncol=nlatent),
            DIFFUSION=matrix(c(diffusion11, 0, 0, 0.3), byrow=TRUE, nrow=nlatent, ncol=nlatent),
            MANIFESTVAR=matrix(c(manifestvar11, 0,0, 0,.3,0, 0,0,.5), nrow=nmanifest, ncol=nmanifest),
            TDPREDEFFECT=matrix(c(tdpredeffect11,0,tdpredeffect12,
              .2,0,-.7), nrow=nlatent, ncol=ntdpred),
            CINT=matrix(c(cint1, 2), nrow=nlatent, ncol=1),
            TDPREDMEANS=TDPREDMEANS,
            T0MEANS=matrix(c(-3,2),nrow=nlatent),
            T0VAR=matrix(c(t0var11,0,0,1),nlatent,nlatent),
            TDPREDVAR=diag(.2,tpoints*ntdpred),
            MANIFESTMEANS=matrix(c(-5,manifestmeans2,1), nrow=nmanifest, ncol=1))
          
          dat=ctGenerate(ctmodelobj=genm, n.subjects=nsubjects, burnin=0, dtmean=1.2, 
            logdtsd=1,simultdpredeffect=TRUE,wide=TRUE)
          # plot(dat[1,1:tpoints])
          ctIndplot(datawide = dat,n.subjects = 1,n.manifest = nmanifest,Tpoints = tpoints)
          
          assign(pars[1,pari],pars[1,pari]) #set free parameter
          
          fitm= ctModel(Tpoints=tpoints,
            n.latent=nlatent, n.manifest=nmanifest, n.TDpred=ntdpred,
            LAMBDA=matrix(c(1, 0,0,1,.6,lambda32), nrow=nmanifest, ncol=nlatent),
            DRIFT=matrix(c(drift11, 1, 0, -0.5), byrow=TRUE, nrow=nlatent, ncol=nlatent),
            DIFFUSION=matrix(c(diffusion11, 0, 0, 0.3), byrow=TRUE, nrow=nlatent, ncol=nlatent),
            MANIFESTVAR=matrix(c(manifestvar11, 0,0, 0,.3,0, 0,0,.5), nrow=nmanifest, ncol=nmanifest),
            TDPREDEFFECT=matrix(c(tdpredeffect11,0,tdpredeffect12,
              .2,0,-.7), nrow=nlatent, ncol=ntdpred),
            CINT=matrix(c(cint1, 2), nrow=nlatent, ncol=1),
            # TDPREDMEANS=TDPREDMEANS,
            T0MEANS=matrix(c(-3,2),nrow=nlatent),
            T0VAR=matrix(c(t0var11,0,0,1),nlatent,nlatent),
            # TDPREDVAR=diag(.2,tpoints*ntdpred),
            MANIFESTMEANS=matrix(c(-5,manifestmeans2,1), nrow=nmanifest, ncol=1))
          
          # dat[1,'TD1_T0']=NA
          ta=proc.time()[3]
          fit=ctFit(dat, fitm,verbose=2,nofit=F,retryattempts=1) #,objective='Kalman')
          tb=proc.time()[3]
          fit=ctCI(ctfitobj = fit,confidenceintervals = pars[1,pari])
          tc=proc.time()[3]
          se=c(fit$mxobj$output$estimate- fit$mxobj$output$standardErrors*2,
            fit$mxobj$output$estimate, 
            fit$mxobj$output$estimate+ 2*fit$mxobj$output$standardErrors,NA,fit$mxobj$output$status$code,tb-ta)
          ci=c(fit$mxobj$output$confidenceIntervals,NA,fit$mxobj$output$status$code,tc-ta)
          if(varpar) ci = exp(ci)
          
          
          # expect_equal(as.numeric(summary(fit)$ctparameters[1]),.4,tolerance=.02)
          
          ldat=ctWideToLong(datawide = dat,Tpoints = tpoints,n.manifest = nmanifest,n.TDpred = ntdpred,n.TIpred = 0)
          ldat=ctDeintervalise(ldat)
          
          sm=ctStanModel(ctmodelobj = fitm,type = 'stanct',indvarying=F)
          
          ta=proc.time()[3]
          sfit=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=400,chains=2,esthyper = F,optimize=F)
          tb=proc.time()[3]
          khdi=c(summary(sfit)$popmeans[1,c('2.5%','50%','97.5%','n_eff','Rhat')],tb-ta)
          
          tmpdir=tempdir()
          tmpdir=gsub('\\','/',tmpdir,fixed=T)
          system(paste0("rm ",tmpdir,'/*.*'))
          
          ta=proc.time()[3]
          sfit=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=400,chains=2,esthyper = F,optimize=F,kalman=F)
          tb=proc.time()[3]
          shdi=c(summary(sfit)$popmeans[1,c('2.5%','50%','97.5%','n_eff','Rhat')],tb-ta)
          
          tmpdir=tempdir()
          tmpdir=gsub('\\','/',tmpdir,fixed=T)
          system(paste0("rm ",tmpdir,'/*.*'))
          
          row = cbind(as.numeric(pars[2,pari]),round(rbind(se, ci,khdi,shdi),3))
          rownames(row)[1]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','stdError')
          rownames(row)[2]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','ci')
          rownames(row)[3]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','kalman hdi')
          rownames(row)[4]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','states hdi')
          
          # expect_equal(as.numeric(ci[2]),as.numeric(hdi[2]),tolerance=.03)
          
          out<-rbind(out,row)
          colnames(out)=c('true','lower','estimate','upper','n_eff','Code/Rhat','seconds')
          print(out)
        }
      }
    }
  }
# }
# )
