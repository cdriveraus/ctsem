if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  
  cores = 2
  
  context('misc')
  test_that("loo", {    
    # for(cores in c(1,cores)){
    #cores=10
    data(AnomAuth)
    
    AnomAuthmodel<-ctModel(LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
      n.latent=2,n.manifest=2, 
      # TIpredNames = 'Y1',
      # DRIFT=c('dr1||||Y1','dr12||||Y1','dr21||||Y1','dr22||||Y1'),
      MANIFESTVAR=diag(0,2),
      Tpoints=5)
    
    aa=ctDeintervalise(ctWideToLong(AnomAuth[1:500,],
      Tpoints = AnomAuthmodel$Tpoints,n.manifest = 2))
    aa[4:20,AnomAuthmodel$manifestNames] <- NA
    
    
    sm <- ctStanModel(AnomAuthmodel,tipredDefault = FALSE)
    sm$pars$indvarying<- FALSE
    
    sf=ctStanFit(aa,
      ctstanmodel = sm, optimize=TRUE,verbose=0,savescores = FALSE,cores=cores)
    sfboot <- ctStanFit(aa,
      ctstanmodel = sm, optimize=TRUE,verbose=0,savescores = FALSE,cores=cores,optimcontrol=list(bootstrapUncertainty=T))
    
    sdat <- sf$standata
    sdat$dokalmanrows[sdat$subject==1] <- 0L #remove 1 subject
    smf <- stan_reinitsf(sf$stanmodel,sdat)
    test_isclose(
      sf$stanfit$optimfit$value - rstan::log_prob(smf,sf$stanfit$rawest), #check ll equiv
      sum(sf$stanfit$transformedparsfull$llrow[sdat$subject==1]))
    
    loo=ctLOO(fit = sf,folds = 10,cores=cores,parallelFolds = T,subjectwise = T)
    loo2=ctLOO(fit = sf,folds = 10,cores=cores,parallelFolds = F,subjectwise = T)
    loo3=ctLOO(fit = sf,folds = sf$standata$nsubjects,cores=cores,parallelFolds = T,subjectwise = T,casewiseApproximation = T)
    loo4=ctLOO(fit = sfboot,folds = sf$standata$nsubjects,cores=cores,parallelFolds = T,subjectwise = T,casewiseApproximation = T)
    
    if(F){ #check casewise approx against leave one subject out
      looFull=ctLOO(fit = sf,folds = sf$standata$nsubjects,cores=20,parallelFolds = T,subjectwise = T,casewiseApproximation = F)
      plot(looFull$outsampleLogLikRow,loo3$outsampleLogLikRow)
      points(looFull$outsampleLogLikRow,loo4$outsampleLogLikRow,col='red')
      abline(0,1)
      
      sum(looFull$outsampleLogLikRow-loo4$outsampleLogLikRow)
      sum(looFull$outsampleLogLikRow-loo3$outsampleLogLikRow)
      
      plot(sf$stanfit$transformedparsfull$llrow[1,]-loo3$outsampleLogLikRow,sf$stanfit$transformedparsfull$llrow[1,]-loo4$outsampleLogLikRow)
    abline(0,1)
    
    sd(apply(loo3$insampleLogLikRowFolds,2,sum,na.rm=T))
    sd(apply(loo$insampleLogLikRowFolds,2,sum,na.rm=T))
    
    l3=melt(data.table(ll=loo3$insampleLogLikeRowFolds,obs=1:nrow(loo3$insampleLogLikeRowFolds)),id.vars='obs')
    l3[,mean:=mean(value,na.rm=T),by=obs]
    l3[,sd:=sd(value,na.rm=T),by=obs]

    
    l=melt(data.table(ll=loo$insampleLogLikRowFolds,obs=1:nrow(loo$insampleLogLikRowFolds)),id.vars='obs')
    l[,mean:=mean(value,na.rm=T),by=obs]
    l[,sd:=sd(value,na.rm=T),by=obs]
    
    ggplot(l3[variable %in% variable[1],],aes(x=obs,y=mean-l[variable %in% variable[1],]$mean))+
      geom_line()+
      # geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd))+
      # geom_line(data=l[variable %in% variable[1],],aes(x=obs,y=mean),color='red')+
      # geom_errorbar(data=l[variable %in% variable[1],],aes(ymin=mean-sd,ymax=mean+sd),color='red')+
      theme_bw()
    
    cor(loo3$outsampleLogLikRow-sf$stanfit$transformedparsfull$llrow[1,],loo$outsampleLogLikRow-sf$stanfit$transformedparsfull$llrow[1,])
    cor(loo2$outsampleLogLikRow-sf$stanfit$transformedparsfull$llrow[1,],loo$outsampleLogLikRow-sf$stanfit$transformedparsfull$llrow[1,])
    
    pari = 6
    range = range(c(loo$foldpars[pari,],loo3$foldpars[pari,]))
    plot(loo$foldpars[pari,],loo3$foldpars[pari,],xlim=range,ylim=range)
    
    }
    
    test_isclose(
      mean(loo2$outsampleLogLikRow-loo$outsampleLogLikRow),
      0,tol=.01)
    
    test_isclose(
      mean(loo3$outsampleLogLikRow-loo$outsampleLogLikRow),
      0,tol=.01)
    
    test_isclose(
    mean(abs(loo3$outsampleLogLikRow -loo$outsampleLogLikRow)),
      0,tol=.05)
    
    test_isclose(
      mean(loo4$outsampleLogLikRow-loo$outsampleLogLikRow),
      0,tol=.01)
    
    test_isclose(
      mean(abs(loo4$outsampleLogLikRow -loo$outsampleLogLikRow)),
      0,tol=.05)
    
    # test_isclose(
    #   mean(c(apply(loo$insampleLogLikRowFolds,1,sum,na.rm=TRUE)/9-sf$stanfit$transformedparsfull$llrow)),
    #   0,
    #   tol=.1)
    # 
    # test_isclose(
    #   mean(c(apply(loo2$insampleLogLikRowFolds,1,sum,na.rm=TRUE)/9-sf$stanfit$transformedparsfull$llrow)),
    #   0,
    #   tol=.1)
    # 
    # 
    # test_isclose(
    #   mean(c(apply(loo3$insampleLogLikRowFolds,1,sum,na.rm=TRUE)/9-sf$stanfit$transformedparsfull$llrow)),
    #   0,
    #   tol=.1)
    
    
    
  })
  
}
