# ctsemtestdat under its name before version 4, so older scripts still run.
# The data are kept once, in ctsemtestdat.rda.
ctstantestdat <- local({
  e <- new.env()
  load('ctsemtestdat.rda', envir = e)
  e$ctsemtestdat
})
