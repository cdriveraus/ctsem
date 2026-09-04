# Regenerate inst/stan/ctsm.stan and inst/stan/ctsmgen.stan.
#
# Those two files are not hand-written source: they are the output of
# ctStanModelWriter(), captured for one reference model so that rstantools can
# precompile them at install time and an ordinary ctFit() need not compile
# anything. When the writer changes, they have to be regenerated, or the
# precompiled program disagrees with the text written for any model that does
# recompile.
#
#   Rscript dev/regenerate-stan.R            # report the diff, write nothing
#   Rscript dev/regenerate-stan.R --write    # regenerate both files
#   Rscript dev/regenerate-stan.R --write --clean-src   # and delete src/
#
# Or from an R session at the package root:
#
#   source("dev/regenerate-stan.R")
#   ctRegenerateStan()                       # dry run
#   ctRegenerateStan(write = TRUE)
#
# The default is a dry run, because a regeneration that is not wanted is
# expensive to notice: the .stan files are large, mostly machine-written, and a
# diff in them reads as noise. Look at the diff first.
#
# --clean-src deletes src/, which forces rstantools to rebuild every stan
# program on the next install. That is a ~20 minute rebuild, and it is only
# needed when the regenerated text must actually be compiled.
#
# This replaces ctsem:::ctsmupdate(), which shipped inside R/ for no reason:
# it was unexported and unreferenced, so only reachable through :::, its
# non-usecurrentwd path pointed at a hardcoded ~/../sync/CT-SEM/ctsem that no
# longer exists, and it asked for confirmation through readline() so it could
# not run non-interactively.

ctRegenerateStan <- function(write = FALSE, clean_src = FALSE, pkgdir = ".") {

  if (!file.exists(file.path(pkgdir, "DESCRIPTION")) ||
      !dir.exists(file.path(pkgdir, "inst/stan"))) {
    stop("pkgdir must be a ctsem package root (needs DESCRIPTION and inst/stan): ",
      normalizePath(pkgdir, mustWork = FALSE), call. = FALSE)
  }
  stanpath <- file.path(pkgdir, "inst/stan")

  if (!"ctsem" %in% loadedNamespaces()) {
    message("loading ctsem (compile = FALSE) ...")
    devtools::load_all(pkgdir, compile = FALSE, quiet = TRUE)
  }

  # The reference model: single subject sunspots, 2 latent / 1 manifest. Its
  # only job is to be a model whose written program exercises the whole writer.
  # forcerecompile = TRUE is load-bearing and not vestigial -- it is what makes
  # the writer emit analytic jacobians rather than the finite-difference
  # fallback, so the text differs without it.
  sunspots <- datasets::sunspot.year
  sunspots <- sunspots[50:(length(sunspots) - (1988 - 1924))]
  id <- 1
  time <- 1749:1924
  datalong <- cbind(id, time, sunspots)

  model <- ctsem::ctModel(type = 'ct', n.latent = 2, n.manifest = 1,
    manifestNames = 'sunspots',
    latentNames = c('ss_level', 'ss_velocity'),
    LAMBDA = matrix(c(1, 'ma1'), nrow = 1, ncol = 2),
    DRIFT = matrix(c(-.0001, 'a21', 1, 'a22'), nrow = 2, ncol = 2),
    MANIFESTMEANS = matrix(c('m1'), nrow = 1, ncol = 1),
    CINT = matrix(c(0, 0), nrow = 2, ncol = 1),
    MANIFESTVAR = diag(.001, 1),
    T0VAR = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
    DIFFUSION = matrix(c(.0001, 0, 0, "diffusion"), ncol = 2, nrow = 2))

  gen <- list()
  for (nm in c("ctsm", "ctsmgen")) {
    txt <- ctsem::ctFit(datalong, model, fit = FALSE,
      gendata = identical(nm, "ctsmgen"), forcerecompile = TRUE)$stanmodeltext
    # Runs of blank lines collapsed to one. Cosmetic, and the committed files
    # have always been written this way.
    gen[[nm]] <- gsub("\\n+", "\n", txt)

    # stanc translates without compiling: cheap, and it catches a writer change
    # that produces text stan cannot parse before anyone spends 20 minutes on
    # a C++ build.
    r <- try(rstan::stanc(model_code = gen[[nm]], verbose = FALSE), silent = TRUE)
    if (inherits(r, "try-error")) {
      stop("generated ", nm, ".stan does not translate:\n", r[1], call. = FALSE)
    }
    message(nm, ".stan: written and translated, ", nchar(gen[[nm]]), " chars")
  }

  changed <- character(0)
  for (nm in names(gen)) {
    f <- file.path(stanpath, paste0(nm, ".stan"))
    old <- if (file.exists(f)) {
      gsub("\n+", "\n", paste0(readLines(f, warn = FALSE), collapse = "\n"))
    } else ""
    oldl <- strsplit(old, "\n", fixed = TRUE)[[1]]
    newl <- strsplit(gen[[nm]], "\n", fixed = TRUE)[[1]]
    if (identical(trimws(oldl), trimws(newl))) {
      message(nm, ".stan: identical to the committed file")
    } else {
      changed <- c(changed, nm)
      message(nm, ".stan: DIFFERS from the committed file (",
        length(oldl), " -> ", length(newl), " lines)")
    }
  }

  if (!write) {
    if (length(changed)) {
      message("dry run. Re-run with write = TRUE (or --write) to update: ",
        paste(changed, collapse = ", "))
    }
    return(invisible(gen))
  }

  for (nm in names(gen)) {
    f <- file.path(stanpath, paste0(nm, ".stan"))
    if (file.exists(f)) file.rename(f, file.path(stanpath, paste0(nm, ".bak")))
    # Binary connection, because cat() to a text connection translates \n to
    # \r\n on Windows and .gitattributes pins these files to LF in the
    # repository and in every checkout. A CRLF working file still commits as
    # LF, so the damage is confined to disk -- which is exactly where a
    # scripted edit matches in one region of a file and silently misses
    # another.
    con <- file(f, open = "wb")
    writeChar(gen[[nm]], con, eos = NULL)
    close(con)
    message("wrote ", f, " (previous kept as ", nm, ".bak)")
  }

  if (clean_src) {
    unlink(file.path(pkgdir, "src"), recursive = TRUE)
    message("removed src/ -- the next install rebuilds every stan program")
  }

  invisible(gen)
}

# Run only when this file is the script Rscript was given, so that
# source()-ing it from another script defines the function and nothing more.
if (any(grepl("regenerate-stan\\.R$", commandArgs(), perl = TRUE))) {
  .a <- commandArgs(trailingOnly = TRUE)
  ctRegenerateStan(write = "--write" %in% .a, clean_src = "--clean-src" %in% .a)
}
