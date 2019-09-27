#unload packages
packs <- c(names(sessionInfo()$otherPkgs), names(sessionInfo()$loadedOnly))
if(length(packs) > 0){ 
  message('Unloading packages -- if any problems occur, please try this from a fresh R session')
  while(length(packs) > 0){
    newpacks <- c()
    for(packi in 1:length(packs)){
      u=try(lapply(paste0('package:', packs), unloadNamespace))
      if(class(u) %in% 'try-error') newpacks <- c(newpacks,packs[packi])
    }
    packs <- newpacks
    Sys.sleep(.1)
  }
}

#install / load build packages
buildpacks <- c('devtools','pkgbuild','remotes')
for(bi in buildpacks){
  if(!requireNamespace(bi)) install.packages(bi)
}
require(pkgbuild)
try(pkgbuild::check_build_tools())

#create / update makevars if needed
cat('Do you already have a MAKEVARS file configured for rstan usage? If unsure, type N')
mv <- readline('Y / N ?')
while(!mv %in% c('Y','N','y','n')) {
  cat('Do you already have a MAKEVARS file configured for rstan usage? If unsure, type N')
  mv <- readline('Y / N ?')
}
if(mv == 'N' || mv =='n'){ #create makevars
  dotR <- file.path(Sys.getenv("HOME"), ".R")
  if (!file.exists(dotR)) dir.create(dotR)
  M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
  if (!file.exists(M)) file.create(M)
  cat("\nCXX14FLAGS=-O3 -mtune=native",
    if( grepl("^darwin", R.version$os)) "CXX14FLAGS += -arch x86_64 -ftemplate-depth-256" else
      if (.Platform$OS.type == "windows") "CXX11FLAGS=-O3 -mtune=native
CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y" else
  "CXX14FLAGS += -fPIC",
    file = M, sep = "\n", append = TRUE)
}

if(.Platform$OS.type == "windows"){
  if(!suppressMessages(pkgbuild::has_rtools())) message('Waiting for Rtools installation to complete...')
  while(! suppressMessages(pkgbuild::has_rtools())){
    Sys.sleep(.1)
  }
}

#check for new versions of critical packages
old = old.packages()
for(importantpack in c('rstan','OpenMx')){
  if(importantpack %in% old) {
    message('Updating ',importantpack)
    install.packages(importantpack,dependencies = TRUE)
  }
}

#install ctsem from github
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true')
remotes::install_github('cdriveraus/ctsem', upgrade='never',INSTALL_opts = "--no-multiarch", 
  dependencies = c("Depends", "Imports"))
