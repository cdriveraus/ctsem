# Setting the julia backend up is the one part of it a user meets before it
# works, so its failure modes decide whether they ever get a fit out of it.
# Everything here runs offline and installs nothing: the URLs are constructed
# rather than fetched, the archive is a synthetic two-file one, and the consent
# checks assert that a session with nobody to ask installs nothing.

test_that("the download URL is built correctly for every supported platform", {
  expected <- c(
    `windows-x86_64` = "winnt/x64/1.12/julia-1.12.7-win64.zip",
    `linux-x86_64` = "linux/x64/1.12/julia-1.12.7-linux-x86_64.tar.gz",
    `linux-aarch64` = "linux/aarch64/1.12/julia-1.12.7-linux-aarch64.tar.gz",
    `macos-x86_64` = "mac/x64/1.12/julia-1.12.7-mac64.tar.gz",
    `macos-aarch64` = "mac/aarch64/1.12/julia-1.12.7-macaarch64.tar.gz")
  # Pinned to the version .ct_julia_archives records hashes for; if the pin
  # moves, these move with it and the hashes have to be refreshed too.
  skip_if_not(identical(ctsem:::.ct_julia_version, "1.12.7"))
  for (platform in names(expected)) {
    archive <- ctsem:::.ctJuliaArchive(ctsem:::.ct_julia_version, platform)
    expect_identical(archive$url,
      paste0("https://julialang-s3.julialang.org/bin/", expected[[platform]]))
    # Only possible because the version is pinned; see .ct_julia_archives.
    expect_match(archive$sha256, "^[0-9a-f]{64}$")
  }
  expect_null(ctsem:::.ctJuliaArchive(platform = "solaris-sparc"))
})

test_that("a version other than the pinned one is downloaded unverified", {
  archive <- ctsem:::.ctJuliaArchive("1.11.9", "linux-x86_64")
  # The minor-version directory has to track the version, or the URL 404s.
  expect_identical(archive$url,
    "https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.9-linux-x86_64.tar.gz")
  expect_true(is.na(archive$sha256))
  expect_message(expect_true(ctsem:::.ctJuliaChecksumOk(tempfile(), NA_character_)),
    "not verified")
})

test_that("the newest Julia is chosen by version, not by name", {
  # "julia-1.9" sorts after "julia-1.12" as a string, which is how a machine
  # with both ends up launching the older one.
  dirs <- file.path("root", c("julia-1.9.4", "julia-1.12.7", "julia-1.10.5"))
  expect_identical(basename(ctsem:::.ctJuliaNewestDir(dirs)), "julia-1.12.7")
  # juliaup's names carry a build suffix, and it leaves half-finished installs
  # behind next to the real ones.
  expect_identical(basename(ctsem:::.ctJuliaNewestDir(
    file.path("root", c("julia-1.12.5+0.x64.w64.mingw32", "julia-temp-t2zkBC")))),
    "julia-1.12.5+0.x64.w64.mingw32")
  expect_null(ctsem:::.ctJuliaNewestDir(character()))
})

test_that("nothing is installed in a session with nobody to ask", {
  skip_if(interactive())
  withr::local_envvar(CTSEM_JULIA_AGREE = "")
  # The default is refusal, which is what keeps R CMD check -- and CRAN's own
  # runs -- from reaching the network.
  expect_false(ctsem:::.ctJuliaAgreed(NULL, "prompt"))
  expect_true(ctsem:::.ctJuliaAgreed(TRUE, "prompt"))
  expect_false(ctsem:::.ctJuliaAgreed(FALSE, "prompt"))
  expect_null(ctsem:::.ctJuliaInstallJulia(agree = FALSE))
})

test_that("consent can be given in advance by environment variable", {
  withr::local_envvar(CTSEM_JULIA_AGREE = "yes")
  expect_true(ctsem:::.ctJuliaAgreed(NULL, "prompt"))
  withr::local_envvar(CTSEM_JULIA_AGREE = "no")
  expect_false(ctsem:::.ctJuliaAgreed(NULL, "prompt"))
  # An explicit argument outranks it, in both directions.
  expect_true(ctsem:::.ctJuliaAgreed(TRUE, "prompt"))
  withr::local_envvar(CTSEM_JULIA_AGREE = "yes")
  expect_false(ctsem:::.ctJuliaAgreed(FALSE, "prompt"))
})

test_that("consent is declined under a Shiny reactive domain, not gated on interactive()", {
  # interactive() is TRUE inside a Shiny server process (stdout is captured for
  # a log pane), which is exactly why .ctJuliaAgreed() must not gate on it
  # directly: askYesNo() would block a server with no console to read from.
  # .ctProgressConsole() already detects a Shiny reactive domain for the same
  # underlying question, so this asserts .ctJuliaAgreed() reaches that check
  # rather than interactive() -- true regardless of whether *this* test session
  # happens to be interactive.
  skip_if_not_installed("shiny")
  withr::local_envvar(CTSEM_JULIA_AGREE = "")
  requireNamespace("shiny", quietly = TRUE)
  testthat::local_mocked_bindings(
    getDefaultReactiveDomain = function() structure(list(), class = "ShinySession"),
    .package = "shiny")
  expect_false(ctsem:::.ctJuliaAgreed(NULL, "prompt"))
})

test_that("declining leaves an error that says how to proceed", {
  skip_if(interactive())
  skip_if_not_installed("JuliaConnectoR")
  withr::local_envvar(CTSEM_JULIA_AGREE = "no")
  expect_error(ctJuliaInstall(force = TRUE), "needs a Julia installation")
  expect_error(ctJuliaInstall(force = TRUE), "ctJuliaInstall\\(agree = TRUE\\)")
})

test_that("unpacking hoists the archive's wrapper directory and validates the result", {
  root <- withr::local_tempdir()
  # The official archives wrap everything in julia-<version>/, so the unpacked
  # layout has to end up at <target>/bin/julia either way. Written as a tarball
  # because R can write one without an external tool; the zip branch differs
  # only in which extractor is called.
  wrapped <- file.path(root, "src", "julia-9.9.9", "bin")
  dir.create(wrapped, recursive = TRUE)
  writeLines("not really julia", file.path(wrapped, ctsem:::.ctJuliaExeName()))
  tarball <- file.path(root, "wrapped.tar.gz")
  withr::with_dir(file.path(root, "src"),
    utils::tar(tarball, "julia-9.9.9", compression = "gzip", tar = "internal"))

  target <- file.path(root, "install", "julia-9.9.9")
  bin <- ctsem:::.ctJuliaUnpack(tarball, target, "tar.gz")
  expect_identical(normalizePath(bin, winslash = "/"),
    normalizePath(file.path(target, "bin"), winslash = "/"))
  expect_true(ctsem:::.ctJuliaIsBinDir(bin))
  # The staging directory is not left behind for a later session to find.
  expect_false(dir.exists(paste0(target, "-partial")))

  # An archive with no julia executable in it is a failed install, rather than
  # an install that fails later, deep inside Pkg.
  dir.create(file.path(root, "empty", "notjulia"), recursive = TRUE)
  writeLines("x", file.path(root, "empty", "notjulia", "readme.txt"))
  bad <- file.path(root, "empty.tar.gz")
  withr::with_dir(file.path(root, "empty"),
    utils::tar(bad, "notjulia", compression = "gzip", tar = "internal"))
  expect_error(ctsem:::.ctJuliaUnpack(bad, file.path(root, "install", "bad"), "tar.gz"),
    "does not contain")
})

test_that("a ctsem-installed Julia is found again without any environment variable", {
  # This is the step that used to cost a JULIA_BINDIR and a restart: a later
  # session has to rediscover the install on its own. R_user_dir() honours
  # R_USER_DATA_DIR, so the search can be pointed somewhere disposable.
  root <- withr::local_tempdir()
  withr::local_envvar(R_USER_DATA_DIR = root, JULIA_BINDIR = "")
  bin <- file.path(ctsem:::.ctJuliaInstallRoot(), "julia-1.12.7", "bin")
  dir.create(bin, recursive = TRUE)
  writeLines("not really julia", file.path(bin, ctsem:::.ctJuliaExeName()))

  expect_identical(normalizePath(ctsem:::.ctJuliaManagedBin(), winslash = "/"),
    normalizePath(bin, winslash = "/"))
  expect_identical(normalizePath(ctsem:::.ctJuliaBin(), winslash = "/"),
    normalizePath(bin, winslash = "/"))
  # ...and an explicitly configured Julia still outranks it.
  other <- file.path(root, "elsewhere", "bin")
  dir.create(other, recursive = TRUE)
  writeLines("not really julia", file.path(other, ctsem:::.ctJuliaExeName()))
  withr::local_envvar(JULIA_BINDIR = other)
  expect_identical(normalizePath(ctsem:::.ctJuliaBin(), winslash = "/"),
    normalizePath(other, winslash = "/"))
})

test_that("a model still prepares on a machine with no Julia at all", {
  # ctFit(fit=FALSE) is pure R, and CRAN's check machines are exactly the case
  # it has to keep working on: JuliaConnectoR installed, because Suggests are,
  # and no Julia anywhere. Asking for the setup here would turn every
  # preparation test into an error rather than a skip.
  skip_if(interactive())
  withr::local_envvar(CTSEM_JULIA_AGREE = "no")
  bin <- ctsem:::.ctJuliaBin
  assignInNamespace(".ctJuliaBin", function(julia_bin = NULL) NULL, ns = "ctsem")
  on.exit(assignInNamespace(".ctJuliaBin", bin, ns = "ctsem"), add = TRUE)

  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = diag(1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diffusion", 1, 1),
    MANIFESTVAR = matrix("residual", 1, 1), MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1)))
  dat <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2), Y1 = 0)
  prepared <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))
  expect_s3_class(prepared, "ctJuliaModel")

  # ...and asking for the fit itself is the thing that says what is missing.
  expect_error(suppressMessages(ctFit(dat, model, backend = "julia")),
    "Julia was not found")
})

test_that("ctJuliaStatus reports rather than errors, and installs nothing", {
  # It is what a user reaches for when the backend will not start, so it has to
  # work on the machine where nothing else does.
  withr::local_envvar(CTSEM_JULIA_AGREE = "no")
  status <- ctJuliaStatus()
  expect_type(status$available, "logical")
  expect_type(status$connectoR, "logical")
  expect_match(status$engine, "^[0-9a-f]{12}$")
})

test_that("ctJuliaSetup declines to instantiate the engine's Julia dependencies without consent", {
  # Pkg.instantiate() is the network-touching step inside ctJuliaSetup(), and on
  # a machine that already has the environment instantiated -- every other test
  # in this suite -- it is never reached at all, because the environment loads
  # on the first try. To exercise the consent gate this test has to simulate an
  # environment that does *not* yet load, and it does that with the smallest
  # possible mock: every call to ctsem's .ctJuliaEval() -- the one route
  # ctJuliaSetup() evaluates Julia code by -- is real except the one it uses to
  # test readiness, which is made to report "not ready". Consent is explicitly
  # declined (CTSEM_JULIA_AGREE="no"), so the
  # function must stop before ever reaching Pkg.instantiate() -- this test
  # cannot itself trigger a download either way.
  skip_without_julia()
  skip_if(interactive())
  withr::local_envvar(CTSEM_JULIA_AGREE = "no")

  # Worded as Julia words a missing package, because only that failure reaches
  # the gate: any other is an engine that is broken, which installing cannot
  # fix, and ctJuliaSetup() says so instead of asking.
  real_juliaEval <- ctsem:::.ctJuliaEval
  testthat::local_mocked_bindings(
    .ctJuliaEval = function(expr) {
      if (identical(expr, "using ContinuousTimeSEM")) {
        stop("simulated: ArgumentError: Package ContinuousTimeSEM not found in ",
          "current path.")
      }
      real_juliaEval(expr)
    },
    .package = "ctsem")

  expect_error(ctJuliaSetup(), "consent was not given")
  expect_error(ctJuliaSetup(), "Julia package dependencies")
})

test_that("a missing package is told apart from a broken engine, in every wording Julia uses", {
  # Recorded from real failures. The parent-package one is what a first fit on
  # a fresh Julia 1.12 depot reports; it went unrecognised, and the user was
  # told installing would not help when installing was the whole remedy.
  missing <- c(
    'LoadError: failed to find source of parent package: "ChainRulesCore"',
    "ArgumentError: Package ContinuousTimeSEM not found in current path.",
    "ArgumentError: Package Optim [429524aa] not found during precompilation",
    "Package ForwardDiff [f6369f11] is required but does not seem to be installed:")
  for (m in missing) expect_true(ctsem:::.ctJuliaLoadFailIsMissing(m), label = m)
  # A docstring error fails at load too, and no install can fix it.
  broken <- c(
    "LoadError: cannot document the following expression:",
    "LoadError: UndefVarError: `optimcontrol` not defined in `ContinuousTimeSEM`")
  for (m in broken) expect_false(ctsem:::.ctJuliaLoadFailIsMissing(m), label = m)
})

# What ctsem keeps under R_user_dir() is removed once outdated, as CRAN asks,
# and never while it may still be in use.
test_that("engine environments go once unused for a month, the current one never", {
  root <- file.path(withr::local_tempdir(), "julia")
  make <- function(name, days_ago) {
    dir <- file.path(root, name)
    dir.create(dir, recursive = TRUE)
    writeLines("x", file.path(dir, "Project.toml"))
    stamp <- file.path(dir, ".ctsem-last-used")
    file.create(stamp)
    Sys.setFileTime(stamp, Sys.time() - days_ago * 86400)
    dir
  }
  current <- make("engine-current", 90)
  stale <- make("engine-stale", 45)
  recent <- make("engine-recent", 3)
  halfcopied <- make("engine-stale-partial", 2)
  other <- make("not-an-engine", 400)
  ctsem:::.ctJuliaPruneEngines(current)
  expect_true(dir.exists(current))
  expect_false(dir.exists(stale))
  expect_true(dir.exists(recent))
  expect_false(dir.exists(halfcopied))
  expect_true(dir.exists(other))
  # Using one is what keeps it: the stamp is touched on every setup.
  ctsem:::.ctJuliaTouchEngine(recent)
  expect_lt(as.numeric(difftime(Sys.time(), file.mtime(file.path(recent,
    ".ctsem-last-used")), units = "mins")), 5)
})

test_that("older ctsem-installed Julias go once a newer one is in", {
  root <- file.path(withr::local_tempdir(), "julia")
  for (name in c("julia-1.10.5", "julia-1.12.7", "julia-1.13.1", "julia-1.12.7-partial",
    "something-else")) {
    dir.create(file.path(root, name, "bin"), recursive = TRUE)
  }
  ctsem:::.ctJuliaPruneInstalls(file.path(root, "julia-1.12.7"))
  expect_setequal(list.files(root), c("julia-1.12.7", "julia-1.13.1", "something-else"))
})
