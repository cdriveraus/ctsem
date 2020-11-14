#' Analyse divergences in a stanfit object
#'
#' @param sf stanfit object.
#' @param nupars either the string 'all', or an integer reflecting how many pars 
#' (from first to nupars) to use.
#'
#' @return A list of four matrices. $locationsort and $sdsort contian the bivariate interactions of 
#' unconstrained parameters, sorted by either the relative location of any divergences, or the relative standard deviation.
#' $locationmeans and $sdmeans collapse across the bivariate interactions to return the means for each parameter.
#' @export
#'
#' @examples
#' \donttest{
#' sunspots<-sunspot.year
#' sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#' id <- 1
#' time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#'
#' #setup model
#' ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
#'  manifestNames='sunspots', 
#'  latentNames=c('ss_level', 'ss_velocity'),
#'  LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
#'  DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
#'  MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
#'  MANIFESTVAR=diag(0,1), #As per original spec
#'  CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'  DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#'
#' #fit
#' ssfit <- ctStanFit(datalong, ssmodel, iter=2, 
#'   optimize=FALSE, chains=1)
#'   
#' stan_checkdivergences(ssfit$stanfit$stanfit) #stan object
#' }
stan_checkdivergences <- function(sf,nupars = 'all'){

samplerps <- get_sampler_params(sf)
if('all' %in% nupars) nupars <- get_num_upars(sf)
e <- rstan::extract(sf)
ea <- as.array(sf) #[,,1:nupars]

ucsnames <- dimnames(ea)$parameters[1:nupars]

ucs2 <- aperm(stan_unconstrainsamples(sf),c(2,1))

#get divergences
dvg <- unlist(lapply(samplerps, function(x) x[(nrow(x)-dim(ea)[1]+1):nrow(x),'divergent__'])) 

#seperate samples into diverging and not
ucsbad <- ucs2[dvg==1,,drop=FALSE]
ucsgood <- ucs2[dvg==0,,drop=FALSE]

#calculate bivariate summary stats
ints <- matrix(NA,nrow=6, ncol= (ncol(ucs2)^2-ncol(ucs2))/2)
count <- 0
cnames <- c()
for(i in 1:(ncol(ucs2)-1)){
  for(j in (i+1):ncol(ucs2)){
    count = count +1
    ibadgc <- (ucsbad[,i]-mean(ucsgood[,i])) * (ucsbad[,j]-mean(ucsgood[,j])) #divergent interactions centered around good interactions
    igoodgc <- (ucsgood[,i]-mean(ucsgood[,i])) * (ucsgood[,j]-mean(ucsgood[,j]))#centered regular interactions centered around good interactions
    ibadbc <- (ucsbad[,i]-mean(ucsbad[,i])) * (ucsbad[,j]-mean(ucsbad[,j])) #divergent interactions centered around bad interactions
    ints[1,count] <- mean(ibadgc)
    ints[2,count] <- mean(igoodgc)
    
    ints[3,count] <- sd(ibadbc)
    ints[4,count] <- sd(igoodgc)
    cnames <- rbind(cnames,c(ucsnames[i], ucsnames[j]))
  }
}

colnames(cnames) <- c('par1','par2')
rownames(ints) <- c('mean_dvg','mean_good','sd_dvg_badc','sd_good_goodc','dvg_rel_location','dvg_rel_sd')

ints[5,] <-  (ints[1,] - ints[2,])/ ints[4,] #location of mean of divergences relative to mean / sd of good samples
ints[6,] <- ints[3,] / ints[4,] #ratio of sd of divergences to sd of good samples

#get ordering of par interactions by stat
badmean <- order(abs(ints[5,]),decreasing = TRUE)
badsd <- order(ints[6,])

ints <- cbind(cnames,t(ints))
sdsort=ints[badsd,c('par1','par2','dvg_rel_location','dvg_rel_sd')]
locationsort = ints[badmean,c('par1','par2','dvg_rel_location','dvg_rel_sd')]

# sdscores <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) sum(c(which(sdsort[,'par1'] %in% x),which(sdsort[,'par2'] %in% x))))))
# locationscores <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) sum(c(which(sdsort[,'par1'] %in% x),which(locationsort[,'par2'] %in% x))))))

sdmeans <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) {
  mean(as.numeric(c(ints[which(ints[,'par1'] %in% x),'dvg_rel_sd'], 
    ints[which(ints[,'par2'] %in% x),'dvg_rel_sd'])))
})))

locationmeans <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) {
  mean(as.numeric(c(ints[which(ints[,'par1'] %in% x),'dvg_rel_location'], 
    ints[which(ints[,'par2'] %in% x),'dvg_rel_location'])))
})))


 sdmeans <- sdmeans[order(sdmeans[,2]),]
 locationmeans <- locationmeans[order(abs(as.numeric(locationmeans[,2])),decreasing = TRUE),]
 
 rownames(sdmeans) <- sdmeans[,1]
 rownames(locationmeans) <- locationmeans[,1]
 sdmeans <- sdmeans[,2,drop=FALSE]
 locationmeans <- locationmeans[,2,drop=FALSE]
 sdmeans[,1] <- as.numeric(sdmeans[,1])
 locationmeans[,1] <- as.numeric(locationmeans[,1])

#output
out<-list(sdmeans=sdmeans,locationsort=locationsort,sdsort=sdsort,locationmeans=locationmeans)

return(out)
}
