skip_on_cran()
skip_on_32bit()
{  # body of the guard this replaced; indentation unchanged


context("knownFits")

#anomauth
test_that("anomauth", {
  
  
  if( .Machine$sizeof.pointer != 4){
    library(ctsem)
    library(testthat)
    
    cores=2
  #library(ctsem);cores=12

  data(AnomAuth)
  AnomAuthmodel<-ctModel(LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
    n.latent=2,n.manifest=2, 
    MANIFESTVAR=diag(0,2),
    Tpoints=5)

   sm1 <- AnomAuthmodel
  sm1$pars$indvarying<- FALSE
  a=Sys.time()
  # sink('bad.txt')
  sf=ctFit(ctDeintervalise(ctWideToLong(AnomAuth,Tpoints = AnomAuthmodel$Tpoints,n.manifest = 2)),
    model= sm1, optimize=TRUE,verbose=0,savescores = FALSE,cores=cores)
  # sink()
  print(Sys.time()-a)
  # PROVENANCE: 23415.929 is the OpenMx -2LL for this model, i.e. what the
  # previous engine produced -- it is `expect_equal(23415.929, AnomAuthfit$
  # mxobj$output$Minus2LogLikelihood)` in ctsem 3cd210cc (Nov 2016) and still
  # is in ctsemOMX's own tests/testthat/test-knownFits.R. It came into the
  # stan tests unchanged, so it is a cross-implementation reference and not a
  # recording of what stan happens to do. The julia backend reproduces it to
  # 5e-5 absolute (test-julia-backend.R).
  test_isclose(23415.929,-2*sf$stanfit$optimfit$value,tol=.01)
  # A stan-side convergence property: `ginfn` is the infinity norm of the
  # gradient at the reported optimum, recorded by the optimiser and, until
  # review J15/R8, asserted nowhere. Its termination *reason* is not the same
  # claim -- both stan optimisers report "no step found" as a termination --
  # so the gradient is what gets checked. See test-stan-convergence.R.
  # Measured 0.016 and 0.023 on two runs of this fit; a fit stalled at its
  # starting values would be orders of magnitude above the bound.
  expect_lt(sf$stanfit$optimfit$ginfn, 1)
  anoms=summary(sf)
  # PROVENANCE: origin unknown -- treat as a regression pin, not ground truth.
  # It entered as `.038 ... tolerance=.004` in a072d351 (Dec 2019) with the
  # stan port, has no counterpart in ctsemOMX (whose summary reports no such
  # column) and no source named anywhere. It pins the reported standard error
  # of the Y1 manifest mean at whatever ctsem produced then.
  test_isclose(.036,anoms$popmeans['mm_Y1','sd'],tol=.01)
 }

})



test_that("oscillator", {
data("Oscillating")

inits <- c(-39.5, -.5, .1, 1, 0, 1, 0.05, .9)
names(inits) <- c("crosseffect","autoeffect", "diffusion",
  "T0var11", "T0var21", "T0var22","m1", "m2")

oscillatingm <- ctModel(type='omx', n.latent = 2, n.manifest = 1, Tpoints = 11, 
  MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1),
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1), 
  T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
  DRIFT = matrix(c(1e-5, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2), 
  CINT = matrix(0, ncol = 1, nrow = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2))#,

  oscillatingm$DRIFT[2,1]="crosseffect|-log1p(exp(-param))-1e-5"
 sm <- ctModelConvertOMX(oscillatingm)
  sm$pars$indvarying<- FALSE
  sf=ctFit(ctDeintervalise(ctWideToLong(Oscillating,Tpoints = oscillatingm$Tpoints,n.manifest = 1)),
    cores=2,verbose=0,
    # optimcontrol=list(carefulfit=T),
    model= sm, optimize=TRUE,savescores = FALSE,priors=FALSE)
  # PROVENANCE: -3461.936 is the OpenMx -2LL for the damped-oscillator example,
  # asserted against `oscillatingf$mxobj$output$Minus2LogLikelihood` in ctsem
  # 3cd210cc (Nov 2016) and still carried, commented out, in ctsemOMX's
  # tests/testthat/test-knownFits.R. Like the AnomAuth value it is what the
  # previous engine produced, so it is a cross-implementation reference rather
  # than a stan recording.
  expect_equal(-3461.936,-2*sf$stanfit$optimfit$value,tolerance=.01)
  # As above: the gradient, not the termination reason. Measured 0.015 and
  # 0.027 on two runs of this fit.
  expect_lt(sf$stanfit$optimfit$ginfn, 1)


})

}
