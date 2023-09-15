if(F){
library(ctsem)
gm <- ctModel(LAMBDA=diag(2), #diagonal factor loading, 2 latents 2 observables
  Tpoints = 5,
  DRIFT=matrix(c(-1,.5,0,-1),2,2), #temporal dynamics
  TRAITVAR = diag(.5,2), #stable latent intercept variance (cholesky factor)
  DIFFUSION=diag(2)) #within person covariance 

ctModelLatex(gm) #to view latex system equations

#when generating data, free pars are set to 0
d <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 100,
  burnin = 20,dtmean = 1))


d$Y2 <- d$Y2 + rnorm(nrow(d),0,.2) #gaussian measurement error
# d$Y2binary <-rbinom(n = nrow(d),size = 1, #create binary data based on the latent
  # prob = ctsem::inv_logit(d$Y2))
d$Y1 <- d$Y1 + rnorm(nrow(d),0,.2) #gaussian measurement error

m <- ctModel(LAMBDA=diag(2),type='stanct',
  manifestNames = c('Y1','Y2'))


f <- ctStanFit(datalong = d,ctstanmodel = m,cores=6,priors = TRUE)
s=summary(f)
s

f2 <- ctStanFit(datalong = d,ctstanmodel = m,cores=6,optimize=F,chains=4,intoverpop = F)
s2=summary(f2)
f3 <- ctStanFit(datalong = d,ctstanmodel = m,cores=6,optimize=F,chains=4,intoverpop = F,intoverstates = F)
s3=summary(f3)
}
