## The feature has only been exercised at k = 2 and k = 6. Does the classifier
## plus rewrite hold at k = 20 (10 DRIFT + 10 DIFFUSION effects on 10 latents)?
## Fit free, so contention-proof: ctIdentify on the specification.
## Enough subjects that the score matrix is not rank limited -- full rank has
## 240 parameters, so 300 subjects.
Sys.setenv(NOT_CRAN='true'); devtools::load_all('.', compile=FALSE, quiet=TRUE)
say <- function(...) { cat('@@', ..., '\n', sep=''); flush.console() }
say('LOADED')
nl <- 10L; k <- 20L; NSUB <- 300L; NT <- 8L
set.seed(11)
Rho <- matrix(0.3, k, k); diag(Rho) <- 1
Z <- matrix(rnorm(NSUB*k), NSUB, k) %*% chol(diag(rep(.5,k)) %*% Rho %*% diag(rep(.5,k)))
o <- vector('list', NSUB)
for(i in seq_len(NSUB)){
  a <- -log1p(exp(Z[i,1:nl])); q <- log1p(exp(Z[i,(nl+1):k]))
  eta <- matrix(0,NT,nl); eta[1,] <- rnorm(nl)
  e <- exp(a); qd <- sqrt(q^2/(-2*a)*(1-e^2))
  for(t in 2:NT) eta[t,] <- e*eta[t-1,] + rnorm(nl,0,qd)
  o[[i]] <- data.frame(id=i, time=seq_len(NT)-1,
    eta + matrix(rnorm(NT*nl,0,sqrt(.2)),NT,nl))
}
dd <- do.call(rbind,o); names(dd)[3:(2+nl)] <- paste0('Y',1:nl)
say('data ', nrow(dd), ' rows, ', NSUB, ' subjects, nlatent ', nl, ', k ', k)
say('full rank population parameters k(k+1)/2 = ', k*(k+1)/2)
say('claimed identified m(m+1)/2 + v*m = ', nl*(nl+1)/2 + nl*nl)
say('claimed flat v(v+1)/2 = ', nl*(nl+1)/2)
mk <- function(){
  DR <- matrix(0,nl,nl); diag(DR) <- paste0('dr',1:nl,'|-log1p_exp(param)')
  DF <- matrix(0,nl,nl); diag(DF) <- paste0('df',1:nl,'|log1p_exp(param)')
  MV <- matrix(0,nl,nl); diag(MV) <- paste0('mv',1:nl,'|log1p_exp(param)')
  m <- suppressMessages(ctModel(type='ct',n.latent=nl,n.manifest=nl,LAMBDA=diag(nl),
    manifestNames=paste0('Y',1:nl), latentNames=paste0('eta',1:nl),
    DRIFT=DR, DIFFUSION=DF, MANIFESTVAR=MV, CINT=matrix(0,nl,1),
    MANIFESTMEANS=matrix(0,nl,1), T0MEANS=matrix(0,nl,1), T0VAR=diag(1,nl)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$matrix %in% c('DRIFT','DIFFUSION') & !is.na(m$pars$param)] <- TRUE
  m
}
## the classifier at scale
roles <- ctsem:::.ctPopEffectRoles(ctsem:::ctModelStatesAndPARS(
  ctsem:::ctModel0DRIFT(mk(), TRUE)$pars, statenames=paste0('eta',1:nl), tdprednames=NULL))
say('classifier: ', nrow(roles), ' effects, ', sum(roles$mean), ' mean-affecting')
## and the prepared specs
for(pr in list(NA,'auto')){
  lab <- if(is.character(pr)) 'auto' else 'NA'
  sp <- try(suppressWarnings(suppressMessages(ctFit(datalong=dd, model=mk(),
    backend='julia', fit=FALSE, intoverpop='augmented', poprank=pr, cores=1L))), silent=TRUE)
  if(inherits(sp,'try-error')){ say('poprank=',lab,' PREPARE FAILED: ',
    substr(conditionMessage(attr(sp,'condition')),1,140)); next }
  pt <- as.data.frame(sp$parameter_table)
  say(sprintf('poprank=%-4s npar %3d  augmented dim %3d', lab,
    ctsem:::.ctBackendNpar(sp), max(pt$row[pt$matrix=='T0MEANS'])))
  id <- try(suppressWarnings(suppressMessages(ctIdentify(dd, mk(), nstart=2L,
    cores=1L, verbose=0L, poprank=pr))), silent=TRUE)
  if(inherits(id,'try-error')) say('   ctIdentify: ',
    substr(conditionMessage(attr(id,'condition')),1,120)) else
    say(sprintf('   ctIdentify npar %3d  nweak %d..%d  rankLimited %s  smallest rel eig %.3g',
      id$npar, id$nweak, id$nweakmax, id$rankLimited, id$smallest))
}
say('DONE')
