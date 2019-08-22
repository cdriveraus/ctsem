if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
library(devtools)
  context("allexamples")
  
  test_that("allexamples1", {
    if (file.exists('../../ctsem/')) run_examples(pkg = '../../ctsem/.') else run_examples(pkg = '../../.')
  })
}
