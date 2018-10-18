stan_reinitsf <- function(model, data){
  suppressMessages(suppressWarnings(suppressOutput(sf<-sampling(model,iter=1,chains=0,init=0,data=data,check_data=FALSE, 
      control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
    return(sf)
}
