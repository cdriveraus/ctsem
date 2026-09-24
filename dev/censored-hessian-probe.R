# Where is ctsem_hessian non-finite on the censored test model?
#   Rscript dev/censored-hessian-probe.R <tree> <lib>
args <- commandArgs(TRUE)
tree <- path.expand(args[1]); lib <- path.expand(args[2])
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(lib, .libPaths()))
suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
ex <- parse(file.path(tree, "tests/testthat/test-julia-censored.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
  startsWith(as.character(e[[2]]), ".")) eval(e, globalenv())
d <- .censored_data(nsubjects = 50, nobs = 10)
m <- .censored_model()
spec <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
  intoverpop = "augmented", fit = FALSE)))
pt <- spec$parameter_table
npar <- max(pt$parnumber, na.rm = TRUE)
names_by <- vapply(seq_len(npar), function(i)
  paste(unique(pt$param[!is.na(pt$parnumber) & pt$parnumber == i]), collapse = "/"),
  character(1))
module <- .ctJuliaModule(spec$project)
obj <- .ctJuliaObjective(structure(spec, class = c("ctJuliaModel", "ctFitModel")))
H_of <- function(fn, x) {
  out <- try(.ctBackendJuliaValue(module[[fn]](obj, .ctJuliaVector(x))), silent = TRUE)
  if (inherits(out, "try-error")) return(structure(NA, err = as.character(out)))
  matrix(as.numeric(out), npar, npar)
}
grad <- function(x) as.numeric(ctJuliaEvaluate(spec, x, gradient = TRUE)$gradient)
fdH <- function(x, h = 1e-5) {
  sapply(seq_len(npar), function(j) {
    up <- x; up[j] <- up[j] + h; dn <- x; dn[j] <- dn[j] - h
    (grad(up) - grad(dn)) / (2 * h)
  })
}
report <- function(label, x) {
  cat("====", label, "\n")
  cat("value", ctJuliaEvaluate(spec, x, gradient = FALSE)$value,
    " gradient finite:", all(is.finite(grad(x))), "\n")
  for (fn in c("ctsem_hessian", "ctsem_hessian_forward")) {
    H <- H_of(fn, x)
    if (length(H) == 1L) { cat(fn, "ERROR:", attr(H, "err"), "\n"); next }
    bad <- which(!is.finite(H), arr.ind = TRUE)
    cat(fn, ": non-finite entries", nrow(bad), "of", length(H), "\n")
    if (nrow(bad)) {
      rows <- sort(unique(c(bad[, 1], bad[, 2])))
      cat("   parameters involved:", paste(rows, names_by[rows], sep = ":", collapse = ", "), "\n")
    }
    assign(paste0("H_", fn), H, envir = globalenv())
  }
  F <- fdH(x)
  cat("fd Hessian finite:", all(is.finite(F)), "\n")
  for (fn in c("ctsem_hessian", "ctsem_hessian_forward")) {
    H <- get0(paste0("H_", fn), envir = globalenv())
    if (is.null(H) || length(H) == 1L) next
    ok <- is.finite(H) & is.finite(F)
    cat(fn, "max |H - fd| over finite entries:", signif(max(abs(H[ok] - F[ok])), 3),
      " scale", signif(max(abs(F[ok])), 3), "\n")
  }
}
set.seed(2)
report("random point (as the gradient test uses)", stats::rnorm(npar, 0, 0.25))
fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
  intoverpop = "augmented", optimcontrol = list(estonly = TRUE, newton = FALSE))))
report("fitted estimate", fit$estimate$raw)
cat("DONE\n")
