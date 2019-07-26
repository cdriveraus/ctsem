if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){

library(ctsem)
  library(devtools)
  run_examples(pkg = '../.')
}
