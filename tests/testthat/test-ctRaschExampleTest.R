if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  context("ctRasch") #develop some expectations here!
  
  test_that("ctRasch1", {
    set.seed( 1234 )
    
    #install software
    # source(file = 'https://github.com/cdriveraus/ctsem/raw/master/installctsem.R')
    
    invlog=function (x) exp(x)/(1 + exp(x))
    n.manifest=7
    
    #gen data
    gm <- ctModel(DRIFT=-.3, DIFFUSION=.3, CINT=.1,
      TRAITVAR=diag(.3,1), #old approach to allow individual variation 
      LAMBDA= rep(1,each=n.manifest),
      n.latent=1,n.manifest=n.manifest,Tpoints=20,
      MANIFESTMEANS=c(0,rep(c(.5,-.5),each=(n.manifest-1)/2)),T0MEANS=-.3,T0VAR=.5)
    
    d=ctGenerate(gm,n.subjects = 20,logdtsd=.2)
    d[,gm$manifestNames] <- rbinom(nrow(d)*gm$n.manifest,size=1,prob=invlog(d[,gm$manifestNames]))
    
    #model to fit
    m <- ctModel( n.latent = 1,
      n.manifest = n.manifest,
      MANIFESTMEANS = c(0,paste0('m',2:n.manifest,'|param|FALSE')), #set prior to N(0,1), disable individual variation
      LAMBDA = rep(1,n.manifest),
      DIFFUSION='diff|log1p_exp(2*param)',
      T0MEANS='t0m|param|TRUE|.2',
      CINT = 'b|param|TRUE|1', #use standard normal for mean prior, individual variation = TRUE (default), default scale for sd
      type = "stanct" )
    
        md <- ctModel( n.latent = 1,
      n.manifest = n.manifest,
      MANIFESTMEANS = c(0,paste0('m',2:n.manifest,'|param|FALSE')), #set prior to N(0,1), disable individual variation
      LAMBDA = rep(1,n.manifest),
      DIFFUSION='diff|log1p_exp(2*param)',
      T0MEANS='t0m|param|TRUE|.2',
      CINT = 'b|param|TRUE|1', #use standard normal for mean prior, individual variation = TRUE (default), default scale for sd
      type = "standt" )
    
    #plot(m)
    
    m$manifesttype[]=md$manifesttype[]=1 #set type to binary
    cores=2
    
        #fit with integration (linearised approximation)
    ro <- ctStanFit( datalong = d,
      ctstanmodel = m,cores=cores,
      # plot=10,verbose=0,
      intoverstates = T,priors=F,
      optimize=T,intoverpop=T)#,optimcontrol=list(stochastic=F))
    so=summary(ro)
    
    
        rod <- ctStanFit( datalong = d,
      ctstanmodel = md,cores=cores,
      # plot=10,verbose=0,
      intoverstates = T,priors=F,
      optimize=T,intoverpop=T)#,optimcontrol=list(stochastic=F))
    sod=summary(rod)
    
    #fit without integration
    r <- ctStanFit( datalong = d,
      #fit=FALSE, #set this to skip fitting and just get the standata and stanmodel objects
      ctstanmodel = m,
      iter = 300,verbose=0,
      control=list(max_treedepth=8),
      priors=TRUE,
      chains = cores,plot=FALSE,
      intoverstates = FALSE,
      optimize=FALSE,intoverpop=FALSE)
    s=summary(r)
    # s
    

    # so
    
    a=cbind(s$popmeans[order(rownames(s$popmeans)),1,drop=FALSE],so$popmeans[order(rownames(so$popmeans)),1])
    colnames(a)=NULL
    # print(a)
    
    testthat::expect_equivalent(
      s$popmeans[order(rownames(s$popmeans)),1],
      so$popmeans[order(rownames(so$popmeans)),1],tol=.1)
    
  })
  
}
