if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  context("ctRasch") #develop some expectations here!
  
  test_that("ctRasch1", {
    set.seed( 1234 )
    invlog=function (x) exp(x)/(1 + exp(x))
    
    gm <- ctModel(DRIFT=-.5, DIFFUSION=.3, CINT=.1,TRAITVAR=diag(.3,1),LAMBDA= c(1,.8,1.2),
      n.latent=1,n.manifest=3,Tpoints=20,
      MANIFESTMEANS=c(0,.5,-.5),T0MEANS=-.3,T0VAR=.5)
    
    d=ctGenerate(gm,n.subjects = 50,logdtsd=0,wide=FALSE)
    # d[,gm$manifestNames] = d[,gm$manifestNames] + rnorm(nrow(d)*gm$n.manifest)
    d[,gm$manifestNames] <- rbinom(nrow(d)*gm$n.manifest,size=1,prob=invlog(d[,gm$manifestNames]))
    
    m <- ctModel( n.latent = 1,
      n.manifest = 3,
      MANIFESTMEANS = c(0,'m2||FALSE','m3||FALSE'),
      LAMBDA = c(1,.8,1.2),
      CINT = 'b',
      type = "stanct" )
    
    m$manifesttype[]=1
    
    r <- ctStanFit( datalong = d,
      ctstanmodel = m,
      iter = 1000,verbose=0,#control=list(max_treedepth=8),
      chains = 3,
      # intoverstates = FALSE,
      # optimize=FALSE,intoverpop=F,
      stationary = FALSE,plot=T )
    s=summary(r)
    s
    
    # m <- ctModel( n.latent = 1,
    #   n.manifest = 3,
    #   # MANIFESTMEANS = c('inv_logit(0+eta1)','inv_logit(m2+eta1)','inv_logit(m3+eta1)'),
    #   # MANIFESTMEANS = c('exp(0+eta1)/(1+exp(0+eta1))','exp(m2+eta1)/(1+exp(m2+eta1))','exp(m3+eta1)/(1+exp(m3+eta1))'),
    #   PARS=c('m2','m3'),
    #   # LAMBDA = c(1,'l2','l3'),
    #   LAMBDA = 0,
    #   CINT = 'b',
    #   type = "stanct" )
    
    # m$manifesttype[]=1
    
    r <- ctStanFit( datalong = d,
      ctstanmodel = m,
      iter = 300,verbose=0,control=list(max_treedepth=8),
      cores = 1,
      # nlcontrol=list(nlmeasurement=T),
      intoverstates = T,nopriors=T,
      optimize=T,intoverpop=T,#fit=F,
      stationary = FALSE)
    s=summary(r)
    s
    
  })
  
}
