if(Sys.getenv("NOT_CRAN")==TRUE & .Machine$sizeof.pointer != 4){

library(ctsem)
  library(devtools)
  run_examples(pkg = 'ctsem')
}
