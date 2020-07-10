.onLoad <- function(libname, pkgname) {
  if(!(.Platform$OS.type=="windows" && .Platform$r_arch=="i386")){
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) loadModule(m, what = TRUE)
  }
}
