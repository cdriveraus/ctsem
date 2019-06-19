stan_reinitsf <- function(model, data,fast=FALSE){
  if(!fast) suppressMessages(suppressWarnings(suppressOutput(sf<-sampling(model,iter=0,chains=0,init=0,data=data,check_data=FALSE, 
    control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
  
  if(fast) sf <- new(model@mk_cppmodule(model),data,0L,getcxxfun(model@dso))
  return(sf)
}

#based on rstan function
getcxxfun <- function(object) {  
  if (length(object@dso_saved) == 0)
    return(function(...) stop("this function should not be called"))
  if (!is_null_cxxfun(object@.CXXDSOMISC$cxxfun)) 
    return(object@.CXXDSOMISC$cxxfun)
  if (!object@dso_saved) 
    stop("the cxx fun is NULL now and this cxxdso is not saved")
  
  # If the file is still loaded  
  # from the help of function dyn.load 
  #   The function dyn.unload unlinks the DLL.  Note that unloading a
  #   DLL and then re-loading a DLL of the same name may or may not
  #   work: on Solaris it uses the first version loaded.
  f <- sub("\\.[^.]*$", "", basename(object@dso_filename)) 
  f2 <- sub("\\.[^.]*$", "", basename(object@.CXXDSOMISC$dso_last_path)) 
  dlls <- getLoadedDLLs()
  if (f2 %in% names(dlls)) { # still loaded 
    DLL <- dlls[[f2]] 
    fx <- cxxfun_from_dll(object@sig, object@.CXXDSOMISC$cxxfun@code, DLL, check_dll = FALSE) 
    assign('cxxfun', fx, envir = object@.CXXDSOMISC) 
    if (!is.null(object@modulename) && object@modulename != '') 
      assign("module", Module(object@modulename, getDynLib(fx)), envir = object@.CXXDSOMISC)
    return(fx) 
  }
  
  # not loaded  
  if (!identical(object@system, R.version$system)) 
    stop(paste("this cxxdso object was created on system '", object@system, "'", sep = ''))
  fx <- cxxfun_from_dso_bin(object) 
  assign('cxxfun', fx, envir = object@.CXXDSOMISC) 
  if (!is.null(object@modulename) && object@modulename != '') 
    assign("module", Module(object@modulename, getDynLib(fx)), envir = object@.CXXDSOMISC)
  return(fx) 
}
