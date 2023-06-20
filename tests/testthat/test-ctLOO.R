if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  
  cores = 2
  
  context('misc')
  test_that("loo", {    
    # for(cores in c(1,cores)){
    #cores=6
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
      ctstanmodel = sm, optimize=TRUE,verbose=0,savescores = FALSE,cores=cores,
      priors=FALSE,
      optimcontrol=list(finishsamples=500,carefulfit=F))
    
    sdat <- sf$standata
    sdat$dokalmanrows[sdat$subject==1] <- 0L #remove 1 subject
    smf <- stan_reinitsf(sf$stanmodel,sdat)
    testthat::expect_equivalent(
      sf$stanfit$optimfit$value - rstan::log_prob(smf,sf$stanfit$rawest), #check ll equiv
      sum(sf$stanfit$transformedparsfull$llrow[sdat$subject==1]))
    
    loo=ctLOO(fit = sf,folds = 10,cores=cores,parallelFolds = T,subjectwise = T)
    loo2=ctLOO(fit = sf,folds = 10,cores=cores,parallelFolds = F,subjectwise = T)
    
    testthat::expect_equivalent(
      sum(loo2$outsampleLogLikRow-loo$outsampleLogLikRow),
      0,tol=1)
    
    testthat::expect_equivalent(
      sum(loo$outsampleLogLikRow),
      sum(sf$stanfit$transformedparsfull$llrow),tol=.1)
    
    testthat::expect_equivalent(
      sum(loo2$outsampleLogLikRow),
      sum(sf$stanfit$transformedparsfull$llrow),tol=.1)
    
  })
  
}
