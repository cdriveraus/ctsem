if(!is.null(names(sessionInfo()$otherPkgs))){ #unload packages
  message('If any problems occur, please try this from a fresh R session')
  suppressWarnings(invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE)))
}

if(!requireNamespace('pkgbuild')) install.packages('pkgbuild')
require(pkgbuild)
try(pkgbuild::check_build_tools())

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

old = old.packages()
for(importantpack in c('rstan','OpenMx')){
  if(importantpack %in% old) {
    message('Updating ',importantpack)
    install.packages(importantpack,dependencies = TRUE)
  }
}

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true')
remotes::install_github('cdriveraus/ctsem', upgrade='never',INSTALL_opts = "--no-multiarch", 
  dependencies = c("Depends", "Imports"))
