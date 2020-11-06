if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
library(ctsem)
library(testthat)
  cores=2

context("knownFits")

#anomauth
test_that("anomauth", {
  
  if( .Machine$sizeof.pointer != 4){

  #cores=6
  data(AnomAuth)
  AnomAuthmodel<-ctModel(LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
    n.latent=2,n.manifest=2, 
    MANIFESTVAR=diag(0,2),
    Tpoints=5)

   sm <- ctStanModel(AnomAuthmodel)
  sm$pars$indvarying<- FALSE
  a=Sys.time()
  # sink('bad.txt')
  sf=ctStanFit(ctDeintervalise(ctWideToLong(AnomAuth,Tpoints = AnomAuthmodel$Tpoints,n.manifest = 2)),
    ctstanmodel = sm, optimize=TRUE,verbose=0,savescores = FALSE,cores=cores,nopriors=TRUE,#forcerecompile = T,
    optimcontrol=list(finishsamples=500,stochastic=T),plot=10,fit=T)
  # sink()
  print(Sys.time()-a)
  testthat::expect_equal(23415.929,-2*sf$stanfit$optimfit$value,tolerance=.01)
  anoms=summary(sf)
  anoms$popmeans['mm_Y1','sd']
  expect_equivalent(.036,anoms$popmeans['mm_Y1','sd'],tolerance=.008)
 }

})



test_that("oscillator", {
data("Oscillating")

inits <- c(-39.5, -.5, .1, 1, 0, 1, 0.05, .9)
names(inits) <- c("crosseffect","autoeffect", "diffusion",
  "T0var11", "T0var21", "T0var22","m1", "m2")

oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11, 
  MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1),
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1), 
  T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
  DRIFT = matrix(c(1e-5, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2), 
  CINT = matrix(0, ncol = 1, nrow = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2))#,
  # startValues = inits)

if( .Machine$sizeof.pointer != 4){
  oscillatingm$DRIFT[2,1]="crosseffect|-log1p(exp(-param))-1e-5"
 sm <- ctStanModel(oscillatingm)
  sm$pars$indvarying<- FALSE
  sf=ctStanFit(ctDeintervalise(ctWideToLong(Oscillating,Tpoints = oscillatingm$Tpoints,n.manifest = 1)),
    cores=2,verbose=0,
    # optimcontrol=list(carefulfit=T),
    ctstanmodel = sm, optimize=TRUE,savescores = FALSE,nopriors=TRUE)
  expect_equal(-3461.936,-2*sf$stanfit$optimfit$value,tolerance=.01)
  
}

})

}
