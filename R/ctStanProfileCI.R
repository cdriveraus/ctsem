ctStanProfileCI <- function(fit, parnames){
  ll=fit$stanfit$optimfit$value
  cores=1
  # fit$standata$profilelltarget=fit$stanfit$optimfit$value
  np=length(fit$stanfit$rawest)
  smf=stan_reinitsf(fit$stanmodel,fit$standata)
  parsouter <- fit$stanfit$rawest
  highpars = fit$stanfit$transformedpars_old[1:np,'97.5%']
  lowpars = fit$stanfit$transformedpars_old[1:np,'2.5%']
  
  parrows <- match(parnames,fit$setup$matsetup$parname) #find raw parameter corresponding to parnames
  cipars <- paste0('rawpopmeans',fit$setup$matsetup$param[parrows])
  # cipars <- fit$setup$matsetup$param[match(parnames,fit$setup$matsetup$parname)]
  
  optimfit=list()
  
  for(pi in 1:length(cipars)){
    for(upperci in c(TRUE,FALSE)){
      if(upperci) init=highpars else init=lowpars
      
      neglpgf<-function(parm) { #used to maximize lp given a fixed parameter
        # print(parm)
        # browser()
        pars<-parsouter
        pars[-which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])] <- parm
        # if(cores > 1 && evaltime > .1){
        #   
        #   out2 <- parallel::clusterApply(cl, 
        #     split(1:standata$nsubjects,sort(1:standata$nsubjects %% min(standata$nsubjects,cores))), function(subjects) parlp(parm,subjects))
        #   out <- try(sum(unlist(out2)),silent=TRUE)
        #   if(standata$verbose > 0) print(out)
        #   attributes(out)$gradient <- try(apply(sapply(out2,function(x) attributes(x)$gradient,simplify='matrix'),1,sum))
        # } else {
        out<-try(log_prob(smf,upars=pars,adjust_transform=TRUE,gradient=TRUE),silent = TRUE)
        # out[1] <- -abs(fit$stanfit$optimfit$value - (out[1]+3.92))
        attributes(out)$gradient <- attributes(out)$gradient[-which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])]
        # }
        
        if(class(out)=='try-error' || is.nan(out)) {
          out[1]=-Inf
          # gradout <<- rep(NaN,length(parm))
        } 
        # print(-out)
        # print(pars)
        return(-out)
      }
      
      
      mizelpginner=list(
        fg=function(innerpars){
          r=neglpgf(innerpars)
          r=list(fn=r[1], gr= -attributes(r)$gradient)
          return(r)
        },
        fn=neglpgf,
        gr=function(innerpars) -attributes(neglpgf(innerpars))$gradient
      )
      
      mizelpgouter_fn <- function(parx){
        message('parx = ',parx)
        parsouter <<-fit$stanfit$rawest #ml inits
        parsouter[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])] <<- parx #parx replaces ml init and stays fixed for inner loop
        pll=abs( (fit$stanfit$optimfit$value - 1.96) - (-mize(parsouter[-which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])],
          fg=mizelpginner,
          max_iter=99999,
          method="L-BFGS",memory=100,
          line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',
          abs_tol=1e-8,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)$f)) + 
          ifelse(upperci && parx < fit$stanfit$rawest[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])], 100, 0) + 
          ifelse(!upperci && parx > fit$stanfit$rawest[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])], 100, 0)
        # print(pll)
        return(pll)
      }
      
      mizelpgouter=list(
        fn=mizelpgouter_fn,
        gr = function(parx){
          parsouter<<-fit$stanfit$rawest #ml inits
          parsouter[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])] <<- parx + 1e-6#parx replaces ml init and stays fixed for inner loop
          up=mizelpgouter_fn(parsouter[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ])
          parsouter[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi])] <<- parx - 1e-6#parx replaces ml init and stays fixed for inner loop
          down=mizelpgouter_fn(parsouter[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ])
          return( (up-down)/(2*1e-6))
        }
      )
  
      optimfit <- mize(init[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ], 
          fg=mizelpgouter, 
          max_iter=99999,
          method="L-BFGS",memory=100,
          line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',
          abs_tol=1e-8,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)$par
      # browser()
      if(upperci) highpars[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ] <- optimfit
      if(!upperci) lowpars[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ] <- optimfit
    }
  }
  
  lowcipars=lowpars[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ]
  highcipars=highpars[which( rownames(fit$stanfit$transformedpars_old[1:np,]) %in% cipars[pi]) ]
  
  low <- sapply(1:length(lowcipars),function(x) tform(param = lowcipars[x],
    transform = fit$setup$matsetup$transform[parrows[x]],
    multiplier = fit$setup$matvalues$multiplier[parrows[x]],
    meanscale = fit$setup$matvalues$meanscale[parrows[x]],
    offset = fit$setup$matvalues$offset[parrows[x]],
    inneroffset = fit$setup$matvalues$inneroffset[parrows[x]])) #constrain_pars(object = smf, lowpars)
  
  high <- sapply(1:length(highcipars),function(x) tform(param = highcipars[x],
    transform = fit$setup$matsetup$transform[parrows[x]],
    multiplier = fit$setup$matvalues$multiplier[parrows[x]],
    meanscale = fit$setup$matvalues$meanscale[parrows[x]],
    offset = fit$setup$matvalues$offset[parrows[x]],
    inneroffset = fit$setup$matvalues$inneroffset[parrows[x]])) 
  
  out <- data.frame(param=parnames,low=low,high=high)
  rownames(out) <- NULL
  
  return(out)
}
