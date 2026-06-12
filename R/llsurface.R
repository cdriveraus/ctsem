# stan_llsurface<- function(sf, sd=2, pars,nsamples=10000){
# 
#   np= length(pars)
#   
#   for(sdi in 1:10){
#   gridi=t(matrix(rnorm(ceiling(nsamples*np/10),pars,sd*(sdi/10^2)),nrow=np))
#   if(sdi ==1) grid=gridi else grid = rbind(grid,gridi)
#   }
#   est=rep(NA,nrow(grid))
#   colnames(grid) = paste0('par',1:np)#rownames(fit$stanfit$transformedpars_old[1:np,])
# 
#   for(ri in 1:nrow(grid)){
#     est[ri] = log_prob(sf,grid[ri,],adjust_transform=TRUE,gradient=FALSE)
#   }
#   
#   
# 
#   dat=data.table::data.table(iter=1:nrow(grid),ll=est, grid)
#   iter <- ll <- value <- NULL
#   dat=data.table::melt(dat,id.vars=c('iter','ll'))
#   
#   ggplot2::ggplot(data = dat, aes(x=value,y=ll))+ geom_smooth() + facet_wrap(. ~ variable,scales = "free")
#     
#   # legend('topright',fit$setup$popsetup$parname[match(fit$setup$popsetup$param,1:ncol(grid))],text.col=col,bty='n')
#   
#   # of=optimizing(object = fit$stanmodel,data= fit$standata,save_iterations=T, hessian=FALSE, iter=1e6, as_vector=FALSE,draws=0,constrained=FALSE,
#   #     tol_obj=tol, tol_rel_obj=0,init_alpha=.00001, tol_grad=0,tol_rel_grad=0,tol_param=0,history_size=50)
#   # of$par$DRIFT
#   # unconstrain_pars(fit$stanfit$stanfit, of$par)
# }


llsurface<- function(fit, sd=2, pars,nsamples=1000,plot=T){

  np= get_num_upars(fit$stanfit$stanfit)
  
  for(sdi in 1:10){
  gridi=t(matrix(rnorm(ceiling(nsamples*np/10),pars,sd*(sdi/10^2)),nrow=np))
  if(sdi ==1) grid=gridi else grid = rbind(grid,gridi)
  }
  est=rep(NA,nrow(grid))
  colnames(grid) = rownames(fit$stanfit$transformedpars_old[1:np,])

  for(ri in 1:nrow(grid)){
    est[ri] = log_prob(fit$stanfit$stanfit,grid[ri,],adjust_transform=TRUE,gradient=FALSE)
  }
  
  

  dat=data.frame(iter=1:nrow(grid),ll=est, grid)
  iter <- ll <- value <- NULL
  dat=melt(dat,id.vars=c('iter','ll'))
  
  if(plot) return(llsurfacePlot(dat))
  else return(dat)
}


llsurfacePlot <- function(dat) ggplot2::ggplot(data = dat, aes(x=value,y=ll))+ geom_smooth() + facet_wrap(. ~ variable,scales = "free")
