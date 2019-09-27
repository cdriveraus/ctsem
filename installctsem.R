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

while(!pkgbuild::check_rtools()){
  Sys.sleep(.1)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true')
  remotes::install_github('cdriveraus/ctsem', INSTALL_opts = "--no-multiarch", 
    dependencies = c("Depends", "Imports"))
}

