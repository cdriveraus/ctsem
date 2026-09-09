# Regenerate man/ and NAMESPACE from the roxygen comments.
#
#   Rscript dev/document.R
#
# Takes about 7 seconds. Use it after any roxygen edit; the .Rd files and
# NAMESPACE are committed, so an edit that is not regenerated never reaches a
# user. That has happened here: several improved @return and @param blocks sat
# in R/ for days while the shipped help pages showed the old text.
#
# Why not devtools::document(). Its default loader calls pkgbuild::build(),
# which compiles the package from a temporary copy with debug flags. That
# rebuild generates far more object sections than the ordinary one and hits the
# mingw assembler's COFF limit:
#
#   as: stanExports_cov.o: too many sections (89954)
#   Fatal error: stanExports_cov.o: file too big
# (verbatim from when cov.stan was still in the tree; the same limit is hit by
# whichever stanExports object the debug rebuild reaches first.)
#
# So it fails after 83 seconds having documented nothing. Same root cause as
# devtools::load_all(compile = TRUE), which cannot succeed here either.
# Documenting does not need a compiler, so the fix is to hand roxygen2 a loader
# that does not use one.
#
# Why pkgload and not roxygen2::load_source. load_source cannot see S3
# registrations, so it rewrites about seven S3method(print, X) lines in
# NAMESPACE into plain export(print.X) and mangles man/print.ctStanModel.Rd.
# It reports success while doing it. The pkgload loader has the real namespace
# and leaves both alone.

setwd(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")))

roxygen2::roxygenise(
  load_code = function(p) pkgload::load_all(p, compile = FALSE, quiet = TRUE)$env)

cat("\nCheck `git diff NAMESPACE` before committing.",
  "\nAn unexpected export() line means a stray @export tag in R/, not a roxygen bug.\n")
