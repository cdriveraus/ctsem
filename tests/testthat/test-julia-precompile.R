# The engine compiles the filter for a fixed set of model shapes at build time,
# so that a first fit does not pay twenty to seventy seconds of specialisation.
# That only works while the shapes it compiled are the ones ctsem's model writer
# actually emits, and nothing about a mismatch is visible: the package loads,
# the fit runs, and it is merely slow again.
#
# The first attempt at the workload hand-wrote a plausible-looking model and got
# four details wrong at once -- `PARS` first rather than last, `JAx` sharing
# DRIFT's parameter numbers, different templates on `T0VAR` and `MANIFESTMEANS`.
# It cost 100 seconds of build time and saved nothing, and only a stopwatch
# said so. The shapes are generated from ctsem itself now
# (tools/generate-precompile-shapes.R); this is what notices when they go stale.
#
# When it fails, regenerate:  Rscript tools/generate-precompile-shapes.R

.precompile_spec <- function(nlatent, intoverpop) {
  set.seed(4)
  dat <- do.call(rbind, lapply(1:3, function(i) {
    d <- data.frame(id = i, time = 0:3)
    for (m in seq_len(nlatent)) d[[paste0("Y", m)]] <- stats::rnorm(4)
    d
  }))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = paste0("Y", seq_len(nlatent)),
    latentNames = paste0("eta", seq_len(nlatent)), LAMBDA = diag(nlatent))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[match(TRUE, model$pars$matrix == "MANIFESTMEANS")] <- TRUE
  suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    fit = FALSE, intoverpop = intoverpop)))
}

test_that("the precompiled shapes still match what ctsem's model writer emits", {
  skip_without_julia()
  module <- ctsem:::.ctJuliaModule(NULL)
  skip_if(is.null(module$ctsem_shape_is_precompiled),
    "engine predates ctsem_shape_is_precompiled")
  nz <- function(x, empty) { x[is.na(x)] <- empty; x }
  v <- ctsem:::.ctJuliaVector

  for (nlatent in 1:2) for (intoverpop in list(TRUE, "laplace")) {
    table <- as.data.frame(.precompile_spec(nlatent, intoverpop)$parameter_table,
      stringsAsFactors = FALSE)
    matched <- ctsem:::.ctBackendJuliaValue(module$ctsem_shape_is_precompiled(
      v(as.character(table$matrix)), v(as.integer(table$row)), v(as.integer(table$col)),
      v(nz(as.integer(table$parnumber), 0L)), v(nz(as.numeric(table$value), NaN)),
      v(nz(as.character(table$transform), "")),
      v(nz(as.character(table$predicttransform), "")),
      v(nz(as.character(table$updatetransform), "")),
      v(nz(as.character(table$tdtransform), ""))))
    expect_true(matched,
      label = sprintf("nlatent=%d intoverpop=%s is a precompiled shape",
        nlatent, as.character(intoverpop)))
  }
})
