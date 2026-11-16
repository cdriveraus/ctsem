# Build ctsem's website (pkgdown) into docs/, and optionally publish it.
#
#   Rscript dev/build-site.R             build docs/ from the tree as it is
#   Rscript dev/build-site.R --publish   build, then replace the gh-pages branch
#                                        on origin with docs/
#
# Run from the package root. Needs pkgdown and rmarkdown, and pandoc (RStudio's
# copy is found). Nothing in the package is installed, loaded or run: the
# reference pages show their examples as code, and the articles are the
# vignettes as they were last rendered.
#
# Why not pkgdown::build_site(). The vignettes fit models with the Julia engine,
# so they ship pre-rendered as vignettes/<name>.html (dev/render-vignettes-
# parallel.R writes them), and pkgdown would render the .qmd sources beside
# them instead, Julia fits and all. So this builds the site from a copy of the
# package in which each <name>.qmd is replaced by <name>.Rmd: the rendered
# page's body converted back to markdown by pandoc, outputs and figures
# included, for pkgdown to lay out as an article. A vignette reaches the site
# when its .html is rendered again and committed.
#
# The copy also gets NEWS.md with pkgdown's heading layout ("# ctsem 4.0.0"
# where the file has "## date" then "### 4.0.0"); the file itself is unchanged.
#
# Publishing pushes one commit, with every file of docs/, to origin's gh-pages
# branch, replacing what was there: the branch holds only the built site, and
# GitHub Pages serves it at the url in _pkgdown.yml. Publish from a clean tree
# at the commit users install (master); the commit message records which.

args <- commandArgs(trailingOnly = TRUE)
publish <- "--publish" %in% args

root <- normalizePath(".", winslash = "/")
if (!file.exists(file.path(root, "_pkgdown.yml")) ||
    !identical(unname(read.dcf(file.path(root, "DESCRIPTION"), "Package")[1, 1]), "ctsem"))
  stop("Run this from the ctsem package root.")

git <- function(..., dir = root) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(dir), ...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0)
    stop("git ", paste(c(...), collapse = " "), " failed:\n", paste(out, collapse = "\n"))
  out
}

if (!rmarkdown::pandoc_available("2.11")) {
  rstudio <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
    "/usr/lib/rstudio/resources/app/bin/quarto/bin/tools")
  found <- rstudio[file.exists(rstudio)]
  if (length(found)) Sys.setenv(RSTUDIO_PANDOC = found[1])
  rmarkdown::find_pandoc(cache = FALSE)
  if (!rmarkdown::pandoc_available("2.11")) stop("pandoc (2.11 or later) is needed.")
}

if (publish) {
  dirty <- git("status", "--porcelain", "--", "man", "vignettes", "README.md",
    "NEWS.md", "NEW-IN-4.md", "_pkgdown.yml", "DESCRIPTION", "inst/CITATION")
  if (length(dirty))
    stop("Commit these before publishing, so the site matches a commit:\n",
      paste(dirty, collapse = "\n"))
}

stage <- file.path(tempfile("ctsem-site-"), "ctsem")
dir.create(file.path(stage, "vignettes"), recursive = TRUE)
dir.create(file.path(stage, "inst"))
for (f in c("DESCRIPTION", "NAMESPACE", "README.md", "NEW-IN-4.md", "_pkgdown.yml"))
  stopifnot(file.copy(file.path(root, f), stage))
stopifnot(file.copy(file.path(root, "inst", "CITATION"), file.path(stage, "inst")))
stopifnot(file.copy(file.path(root, "man"), stage, recursive = TRUE))
if (dir.exists(file.path(root, "pkgdown")))
  stopifnot(file.copy(file.path(root, "pkgdown"), stage, recursive = TRUE))

# NEWS.md: "## <date>" then "### <version>" becomes "# ctsem <version>" with
# the date beneath it, which is the layout pkgdown reads versions from.
news <- paste(readLines(file.path(root, "NEWS.md"), encoding = "UTF-8"), collapse = "\n")
news <- sub("^# ctsem News\\s*\n", "", news)
nversions <- lengths(regmatches(news, gregexpr("(?m)^###", news, perl = TRUE)))
news <- gsub("(?m)^## +(\\S+)\\s*\n(?:\\s*\n)*### *(?:ctsem +)?(\\S+)[ \t]*$",
  "# ctsem \\2\n\n*\\1*", news, perl = TRUE)
news <- gsub("(?m)^### *ctsem +(\\S+)[ \t]*$", "# ctsem \\1", news, perl = TRUE)
if (lengths(regmatches(news, gregexpr("(?m)^# ctsem ", news, perl = TRUE))) != nversions ||
    grepl("(?m)^##", news, perl = TRUE))
  stop("NEWS.md has a heading this script does not expect; see its NEWS section.")
writeLines(enc2utf8(news), file.path(stage, "NEWS.md"), useBytes = TRUE)

# Pandoc writes a bracket in text as \[, which rmarkdown reads back as display
# math, so an interval such as [0.21, 0.48] became an equation. This filter
# writes brackets in text as entities instead; code is not text and keeps its.
brackets <- tempfile(fileext = ".lua")
writeLines(c(
  "function Str(el)",
  "  if not el.text:find('[%[%]]') then return nil end",
  "  local out, rest = pandoc.List(), el.text",
  "  while #rest > 0 do",
  "    local i = rest:find('[%[%]]')",
  "    if not i then out:insert(pandoc.Str(rest)); break end",
  "    if i > 1 then out:insert(pandoc.Str(rest:sub(1, i - 1))) end",
  "    out:insert(pandoc.RawInline('markdown', rest:sub(i, i) == '[' and '&#91;' or '&#93;'))",
  "    rest = rest:sub(i + 1)",
  "  end",
  "  return out",
  "end"), brackets)

# The articles: every vignette that ships pre-rendered (an .html.asis stub).
stubs <- list.files(file.path(root, "vignettes"), "\\.html\\.asis$")
if (!length(stubs)) stop("No pre-rendered vignettes found in vignettes/.")
for (stub in stubs) {
  name <- sub("\\.html\\.asis$", "", stub)
  page <- file.path(root, "vignettes", paste0(name, ".html"))
  entry <- grep("VignetteIndexEntry", readLines(file.path(root, "vignettes", stub)), value = TRUE)
  title <- sub(".*\\\\VignetteIndexEntry\\{(.*)\\}.*", "\\1", entry)
  if (!file.exists(page) || length(title) != 1) stop("vignettes/", stub, ": no page or no title.")
  html <- paste(readLines(page, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  # pkgdown sets the title and the contents itself.
  html <- sub("(?s)<header id=\"title-block-header\">.*?</header>", "", html, perl = TRUE)
  html <- sub("(?s)<nav id=\"TOC\".*?</nav>", "", html, perl = TRUE)
  # Sections are h2 under pkgdown's h1 title; some vignettes write them as h1.
  if (!grepl("<h[12][ >]", html)) stop(name, ": no sections found in its page.")
  shift <- if (grepl("<h1[ >]", html)) "--shift-heading-level-by=1"
  tmp <- tempfile(fileext = ".html")
  writeLines(enc2utf8(html), tmp, useBytes = TRUE)
  md <- tempfile(fileext = ".md")
  rmarkdown::pandoc_convert(tmp, to = "markdown", from = "html", output = md,
    options = c("--wrap=none", paste0("--lua-filter=", brackets), shift))
  body <- readLines(md, encoding = "UTF-8")
  # knitr reads `r ...` anywhere outside a chunk as R code to run.
  if (any(grepl("`r[ #]", body))) stop(name, ": inline `r` code would be run by knitr.")
  writeLines(enc2utf8(c("---", paste0("title: \"", gsub("\"", "\\\\\"", title), "\""), "---", "", body)),
    file.path(stage, "vignettes", paste0(name, ".Rmd")), useBytes = TRUE)
}

docs <- file.path(root, "docs")
pkg <- pkgdown::as_pkgdown(stage, override = list(destination = docs))
if (dir.exists(docs)) pkgdown::clean_site(pkg, quiet = TRUE)
pkgdown::build_site(pkg, examples = FALSE, new_process = FALSE, install = FALSE,
  preview = FALSE, override = list(destination = docs))
stopifnot(file.copy(file.path(root, "vignettes", "hierarchicalmanual.pdf"),
  file.path(docs, "articles"), overwrite = TRUE))
built <- file.path(docs, "articles", paste0(sub("\\.html\\.asis$", "", stubs), ".html"))
if (!all(file.exists(built))) stop("Articles missing from docs/: ", paste(built[!file.exists(built)], collapse = ", "))
# An article's "Source:" link names the .Rmd made above; the source is the .qmd.
for (page in built) {
  name <- sub("\\.html$", "", basename(page))
  html <- readLines(page, encoding = "UTF-8", warn = FALSE)
  rmd <- paste0("vignettes/", name, ".Rmd")
  if (!any(grepl(rmd, html, fixed = TRUE))) stop(page, ": no source link to ", rmd, ".")
  writeLines(enc2utf8(gsub(rmd, paste0("vignettes/", name, ".qmd"), html, fixed = TRUE)), page, useBytes = TRUE)
}
cat("Built", docs, "with", length(built), "articles.\n")

if (publish) {
  head <- git("rev-parse", "HEAD")
  branch <- git("rev-parse", "--abbrev-ref", "HEAD")
  remote <- git("remote", "get-url", "origin")
  site <- tempfile("ctsem-pages-")
  dir.create(site)
  stopifnot(all(file.copy(list.files(docs, full.names = TRUE, all.files = TRUE, no.. = TRUE),
    site, recursive = TRUE)))
  file.create(file.path(site, ".nojekyll"))
  who <- c("-c", paste0("user.name=", shQuote(git("config", "user.name"))),
    "-c", paste0("user.email=", shQuote(git("config", "user.email"))))
  git("init", "-q", "-b", "gh-pages", dir = site)
  git("add", "-A", dir = site)
  git(who, "commit", "-q", "-m", shQuote(paste0("Site built from ", branch, " ", head)), dir = site)
  git("push", "--force", shQuote(remote), "gh-pages", dir = site)
  cat("Published ", substr(head, 1, 8), " (", branch, ") to gh-pages on ", remote, ".\n", sep = "")
}
