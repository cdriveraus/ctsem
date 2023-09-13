
# stan_confidenceRegion <-function(stanfit,parstrings,prefuncstring='(', joinfuncstring=' + ',postfuncstring=')'){
#   mc=As.mcmc.list(stanfit)
#   mc=do.call(rbind,mc)
#   
#   pars <- lapply(parstrings,function(x) paste0(colnames(mc)[grep(x,colnames(mc),fixed=TRUE)]))
#   parsref <- lapply(parstrings,function(x) paste0('mc[,"',colnames(mc)[grep(x,colnames(mc),fixed=TRUE)],'"]'))
#   
#   # if(length(parstrings) > 1 & matchingindices==TRUE & !all(lapply(pars,length)==length(pars[[1]]))) stop ('matchingindices=TRUE but unequal numbers of parameters found matching parstrings')
#   
#   for(pari in 1:length(pars[[1]])){
#     a <- cbind(eval(parse(text=paste0(prefuncstring, paste0(lapply(parsref,function(x) x[pari]),collapse=joinfuncstring),postfuncstring))))
#     colnames(a) <- paste0(prefuncstring, paste0(lapply(pars,function(x) x[pari]),collapse=joinfuncstring),postfuncstring)
#     if(pari==1) out<-a else out <- cbind(out,a)
#   }
#   
#   return(out)
# }
