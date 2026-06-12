if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  
  context("ctRasch") #develop some expectations here!
  
  test_that("ctRasch1", {
    set.seed( 1234 )
    cores=2

    invlog=function (x) exp(x)/(1 + exp(x))
    n.manifest=7
    
    #gen data
    nsubjects <- 20
    cint <- rnorm(nsubjects, mean = .1, sd = .3)
    gm <- ctModel(DRIFT=-.3, DIFFUSION=.3, CINT=.1,
      LAMBDA= rep(1,each=n.manifest),
      n.latent=1,n.manifest=n.manifest,Tpoints=20,
      MANIFESTMEANS=c(0,rep(c(.5,-.5),each=(n.manifest-1)/2)),T0MEANS=-.3,T0VAR=.5)
    
    dlist <- vector("list", nsubjects)
    for(i in seq_len(nsubjects)){
      gm_i <- gm
      gm_i$CINT[] <- cint[i]
      d_i <- ctGenerate(gm_i, n.subjects = 1, logdtsd = .2)
      d_i[, "id"] <- i
      dlist[[i]] <- d_i
    }
    d <- do.call(rbind, dlist)
    d[,gm$manifestNames] <- rbinom(nrow(d)*gm$n.manifest,size=1,prob=invlog(d[,gm$manifestNames]))
    
    #model to fit
    m <- ctModel( n.latent = 1,
      n.manifest = n.manifest,
      MANIFESTMEANS = c(0,paste0('m',2:n.manifest,'|param|FALSE')), #set prior to N(0,1), disable individual variation
      LAMBDA = rep(1,n.manifest),
      DIFFUSION='diff|log1p_exp(2*param)',
      T0MEANS='t0m|param|TRUE|.2',
      CINT = 'b|param|TRUE|1', #use standard normal for mean prior, individual variation = TRUE (default), default scale for sd
      type = "ct" )
    
        md <- ctModel( n.latent = 1,
      n.manifest = n.manifest,
      MANIFESTMEANS = c(0,paste0('m',2:n.manifest,'|param|FALSE')), #set prior to N(0,1), disable individual variation
      LAMBDA = rep(1,n.manifest),
      DIFFUSION='diff|log1p_exp(2*param)',
      T0MEANS='t0m|param|TRUE|.2',
      CINT = 'b|param|TRUE|1', #use standard normal for mean prior, individual variation = TRUE (default), default scale for sd
      type = "dt" )

    m$manifesttype[]=md$manifesttype[]=1 #set type to binary
    
        #fit with integration (linearised approximation)
    ro <- ctFit( datalong = d,
      ctstanmodel = m,cores=cores,
      # plot=10,verbose=0,
      intoverstates = T,priors=F,
      optimize=T,intoverpop=T)#,optimcontrol=list(stochastic=F))
    so=summary(ro)
    
    
        rod <- ctFit( datalong = d,
      ctstanmodel = md,cores=cores,
      # plot=10,verbose=0,
      intoverstates = T,priors=F,
      optimize=T,intoverpop=T)#,optimcontrol=list(stochastic=F))
    sod=summary(rod)
    
    #fit without integration
    r <- ctFit( datalong = d,
      #fit=FALSE, #set this to skip fitting and just get the standata and stanmodel objects
      ctstanmodel = m,
      iter = 300,verbose=0,
      control=list(max_treedepth=8),
      priors=TRUE,
      chains = cores,plot=FALSE,
      intoverstates = FALSE,
      optimize=FALSE,intoverpop=FALSE)
    s=summary(r)

    a=cbind(s$popmeans[order(rownames(s$popmeans)),1,drop=FALSE],so$popmeans[order(rownames(so$popmeans)),1])
    colnames(a)=NULL
    # print(a)
    
    test_isclose(
      s$popmeans[order(rownames(s$popmeans)),1],
      so$popmeans[order(rownames(so$popmeans)),1],tol=.2)
    
  })
  
}
