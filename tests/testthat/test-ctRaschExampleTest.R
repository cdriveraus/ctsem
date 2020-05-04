if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  context("ctRasch") #develop some expectations here!
  
  test_that("ctRasch1", {
    set.seed( 1234 )
    invlog=function (x) exp(x)/(1 + exp(x))
    
    gm <- ctModel(DRIFT=-.3, DIFFUSION=.3, CINT=.1,TRAITVAR=diag(.3,1),LAMBDA= c(1,.8,1.2),
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
      iter = 100,verbose=0,control=list(max_treedepth=8),
      chains = 2,
      intoverstates = TRUE,
      optimize=FALSE,intoverpop=F,
      stationary = FALSE)
    s=summary(r)
    s
    
    ro <- ctStanFit( datalong = d,
      ctstanmodel = m,cores=1,
      iter = 300,verbose=0,control=list(max_treedepth=8),
      intoverstates = T,nopriors=T,
      optimcontrol = list(stochastic=T),
      optimize=T,intoverpop=T,#fit=F,
      stationary = FALSE)
    so=summary(ro)
    so
    
  })
  
}
