if(F){
library(ctsem)
gm <- ctModel(LAMBDA=diag(2), #diagonal factor loading, 2 latents 2 observables
  Tpoints = 7,
  DRIFT=matrix(c(-1,.5,0,-1),2,2), #temporal dynamics
  MANIFESTVAR=diag(.2,2), #measurement error
  DIFFUSION=diag(2)) #within person covariance 

ctModelLatex(gm) #to view latex system equations

#when generating data, free pars are set to 0
nsubjects <- 100
traitChol <- diag(.5,2)
subjectCint <- t(replicate(nsubjects, as.numeric(traitChol %*% rnorm(2))))
dlist <- vector("list", nsubjects)
for(i in seq_len(nsubjects)){
  gm_i <- gm
  gm_i$CINT <- matrix(subjectCint[i,], ncol = 1)
  d_i <- ctGenerate(ctmodelobj = gm_i,n.subjects = 1,logdtsd = .1,burnin = 20,dtmean = 1)
  d_i[, "id"] <- i
  dlist[[i]] <- d_i
}
d <- data.frame(do.call(rbind, dlist))


d$Y2 <- d$Y2 + rnorm(nrow(d),0,.2) #gaussian measurement error
# d$Y2binary <-rbinom(n = nrow(d),size = 1, #create binary data based on the latent
  # prob = ctsem::inv_logit(d$Y2))
d$Y1 <- d$Y1 + rnorm(nrow(d),0,.2) #gaussian measurement error

m <- ctModel(LAMBDA=diag(2),type='ct',TRAITVAR='auto',Tpoints=7,
  manifestNames = c('Y1','Y2'))
system.time({f <- ctFit(datalong = d,ctstanmodel = m,cores=2,priors = F)})

# library(ctsemOMX)
# system.time({fo<-ctsemOMX::ctFit(ctmodelobj = m,dat = d,stationary = NULL)})

s=summary(f)
s

f2 <- ctFit(datalong = d,ctstanmodel = m,cores=6,optimize=F,chains=4,intoverpop = F)
s2=summary(f2)
f3 <- ctFit(datalong = d,ctstanmodel = m,cores=6,optimize=F,chains=4,intoverpop = F,intoverstates = F)
s3=summary(f3)
}
