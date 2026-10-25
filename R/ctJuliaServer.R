# The R bridge's Julia server, precompiled ------------------------------------
#
# A fresh session's first fit spent about six seconds compiling, and almost
# none of it was the engine's: JuliaConnectoR's `main.jl` includes its server
# (`RConnector`, 2265 lines) into `Main` at run time, so every session compiled
# the server's methods (3.9 s) and Base's for its types (1.9 s), and `main.jl`
# imported Pkg and read the source again (0.9 s); a package image holds only
# packages (review/PRECOMPILE-connector-proposal-2026-10-25.md). So ctsem builds
# the installed JuliaConnectoR's own server files into a small package,
# CtsemRServer (inst/julia/rserver/), compiles it once, and starts its sessions
# from that package; JuliaConnectoR connects to the server it finds at
# `JULIACONNECTOR_SERVER`, as it does to any server, and from there the session
# is the one `main.jl` would have started. Anything that does not go to plan --
# a build that fails, a server that does not answer, JuliaConnectoR internals
# that have moved -- falls back to JuliaConnectoR's own start, quietly: the cost
# is the old start, never a failed session.
#
# The package lives under R_user_dir("ctsem", "cache") in a directory keyed by
# the server's source -- the installed JuliaConnectoR's Julia files and ctsem's
# template -- and the Julia version, so an update to either builds a new one
# and a stale server never runs. Its environment holds the real Tables.jl,
# resolved offline from what JuliaConnectoR installed: JuliaConnectoR checks
# for its dummy Tables on connecting and restarts a session that has it.

# JuliaConnectoR's start internals, or NULL where they have moved.
.ctJuliaServerInternals <- function() {
  tryCatch({
    ns <- asNamespace("JuliaConnectoR")
    out <- list(exe = get("getJuliaExecutablePath", envir = ns, inherits = FALSE),
      env = get("getJuliaEnv", envir = ns, inherits = FALSE),
      files = system.file("Julia", package = "JuliaConnectoR"))
    if (!is.function(out$exe) || !is.function(out$env) ||
        !file.exists(file.path(out$files, "RConnector.jl"))) NULL else out
  }, error = function(e) NULL)
}

# How the package is built from its sources: raised with any change to the
# layout or the build below, so a package built the old way is never taken for
# one built this way.
.ctJuliaServerRecipe <- 2L

# Where the package for this server source, recipe and Julia lives.
.ctJuliaServerDir <- function(internals, bin) {
  template <- system.file("julia", "rserver", package = "ctsem")
  files <- c(sort(list.files(internals$files, full.names = TRUE)),
    sort(list.files(template, full.names = TRUE)))
  if (!nzchar(template) || !length(files)) return(NULL)
  manifest <- c(paste(basename(files), tools::md5sum(files)),
    paste("julia", .ctJuliaCachedVersion(bin)),
    paste("recipe", .ctJuliaServerRecipe), deparse(.ctJuliaServerPatches))
  handle <- tempfile(fileext = ".txt")
  on.exit(unlink(handle), add = TRUE)
  writeLines(manifest, handle)
  key <- substr(unname(tools::md5sum(handle)), 1L, 12L)
  gsub("\\", "/", file.path(tools::R_user_dir("ctsem", which = "cache"), "julia",
    paste0("rserver-", key)), fixed = TRUE)
}

# The built package's directory, building it on first use; NULL where it
# cannot be built. A build is a one-off `julia` process: the package's files
# laid out, then `Pkg.develop` into its own environment, offline, and `using`,
# which compiles its image. A directory without `ready` is a build in progress
# elsewhere or one that failed (`failed`, with Julia's output): either way this
# session starts the old way, and a failed build is not retried until its key
# changes or a week has passed.
.ctJuliaServerReady <- function(internals, bin) {
  dir <- .ctJuliaServerDir(internals, bin)
  if (is.null(dir)) return(NULL)
  if (file.exists(file.path(dir, "ready"))) return(dir)
  if (dir.exists(dir)) {
    failed <- file.path(dir, "failed")
    age <- as.numeric(difftime(Sys.time(), file.mtime(dir), units = "days"))
    building <- !file.exists(failed) && isTRUE(age < 1 / 48)
    if (building || (file.exists(failed) && isTRUE(age < 7))) return(NULL)
    unlink(dir, recursive = TRUE)
  }
  staging <- paste0(dir, "-partial-", Sys.getpid())
  unlink(staging, recursive = TRUE)
  laid <- tryCatch(.ctJuliaServerLayout(staging, internals), error = function(e) FALSE)
  if (!isTRUE(laid) || !isTRUE(suppressWarnings(file.rename(staging, dir)))) {
    unlink(staging, recursive = TRUE)
    return(NULL)
  }
  message("Compiling the R-Julia bridge's server, once for this JuliaConnectoR and Julia ...")
  # Its own width, as an engine build's (`.ctJuliaBuildWidth()`): Julia writes
  # an image on a thread per CPU and Pkg builds that many packages at once.
  width <- as.character(.ctJuliaOr(.ct_julia_cache$build_cores, .ctJuliaDefaultCores()))
  saved <- Sys.getenv(c("JULIA_IMAGE_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
    "OPENBLAS_NUM_THREADS"), unset = NA)
  for (k in names(saved)) if (is.na(saved[[k]])) do.call(Sys.setenv, stats::setNames(list(width), k))
  on.exit(for (k in names(saved)) if (is.na(saved[[k]])) Sys.unsetenv(k), add = TRUE)
  script <- file.path(dir, "build.jl")
  writeLines(c(
    "using Pkg",
    "Pkg.offline(true)",
    sprintf("Pkg.activate(%s; io = devnull)", .ctJuliaString(file.path(dir, "env"))),
    # Relative, so the environment's manifest does not name this machine's path.
    "cd(() -> Pkg.develop(Pkg.PackageSpec(path = \"../CtsemRServer\"); io = devnull),",
    sprintf("  %s)", .ctJuliaString(file.path(dir, "env"))),
    # Loaded from its image, not merely loaded: a package whose image fails
    # to build still loads, uncompiled, and every session would pay for both.
    "using CtsemRServer",
    "Base.isprecompiled(Base.identify_package(\"CtsemRServer\")) &&",
    "  print(\"CTSEM-RSERVER-READY\")"), script)
  out <- tryCatch(suppressWarnings(system2(file.path(bin, .ctJuliaExeName()),
    c("--startup-file=no", "--history-file=no", shQuote(script)),
    stdout = TRUE, stderr = TRUE, timeout = 900)), error = function(e) conditionMessage(e))
  if (any(grepl("CTSEM-RSERVER-READY", out, fixed = TRUE))) {
    file.create(file.path(dir, "ready"))
    .ctJuliaServerPrune(dir)
    return(dir)
  }
  writeLines(out, file.path(dir, "failed"))
  NULL
}

# What the server's source needs to be a package, each as `(file, from, to)`
# with `from` found exactly once, or the layout fails and the session starts the
# old way: a JuliaConnectoR whose source has moved on is never edited blind.
#
# - Its Tables is the package's, not `Main`'s, which a package cannot import.
# - `ElementList(unnamed::Vector{Any})` is defined twice in reading.jl, by two
#   methods' optional arguments, and a package image refuses a method
#   overwritten while it is built. Both give `names` empty, so the second's
#   `named` loses its default and the first's stands, which is what the
#   overwrite left anyway.
.ctJuliaServerPatches <- list(
  list("RConnector.jl", "import Main.Tables: istable", "import ..Tables: istable"),
  list("reading.jl", paste0("function ElementList(unnamed::Vector{Any},\n",
    "      named::Dict{Symbol, Any} = Dict{Symbol, Any}(),"),
    paste0("function ElementList(unnamed::Vector{Any},\n",
      "      named::Dict{Symbol, Any},")))

# The package's files, laid out under `dir`: JuliaConnectoR's server source as
# a package (`.ctJuliaServerPatches`) beside ctsem's template, and the
# server's start script.
.ctJuliaServerLayout <- function(dir, internals) {
  template <- system.file("julia", "rserver", package = "ctsem")
  src <- file.path(dir, "CtsemRServer", "src")
  dir.create(file.path(src, "connector"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "env"), showWarnings = FALSE)
  copied <- file.copy(list.files(internals$files, full.names = TRUE),
    file.path(src, "connector"))
  if (!all(copied)) return(FALSE)
  for (patch in .ctJuliaServerPatches) {
    path <- file.path(src, "connector", patch[[1L]])
    text <- gsub("\r\n", "\n", paste(readLines(path, warn = FALSE), collapse = "\n"),
      fixed = TRUE)
    if (lengths(regmatches(text, gregexpr(patch[[2L]], text, fixed = TRUE))) != 1L) {
      return(FALSE)
    }
    writeLines(sub(patch[[2L]], patch[[3L]], text, fixed = TRUE), path)
  }
  file.copy(file.path(template, c("CtsemRServer.jl", "precompile_statements.jl")), src)
  file.copy(file.path(template, "serve.jl"), dir)
  writeLines(c(
    "name = \"CtsemRServer\"",
    "uuid = \"a2cb8021-53d5-4531-812c-902b21fa8bbe\"",
    "version = \"0.1.0\"",
    "",
    "[deps]",
    "InteractiveUtils = \"b77e0a4c-d291-57a0-90e8-8db25a27a240\"",
    "Logging = \"56ddb016-857b-54e1-b83d-db4d58db5568\"",
    "Pkg = \"44cfe95a-1eb2-52ea-b672-e2afdf69b78f\"",
    "REPL = \"3fa0cd96-eef1-5676-8a61-b3b8758bbffb\"",
    "Sockets = \"6462fe0b-24de-5631-8697-dd941f90decc\"",
    "Tables = \"bd369af6-aec1-5ad0-b16a-f7cc5008161c\""),
    file.path(dir, "CtsemRServer", "Project.toml"))
  TRUE
}

# Other servers' directories, once this one is built: a JuliaConnectoR or
# Julia update leaves the old one unused. A directory another R session's
# Julia has open may refuse on Windows, and is left.
.ctJuliaServerPrune <- function(keep) {
  dirs <- list.dirs(dirname(keep), recursive = FALSE)
  old <- dirs[grepl("^rserver-", basename(dirs)) & basename(dirs) != basename(keep)]
  for (d in old) try(unlink(d, recursive = TRUE), silent = TRUE)
  invisible(NULL)
}

# Start a session from the precompiled server and connect JuliaConnectoR to it,
# as JuliaConnectoR's `runJuliaServer()` and `juliaConnection()` would with
# `main.jl`: the same options (`JULIACONNECTOR_JULIAOPTS`), the same port file
# and single client, the same Julia and environment. TRUE once connected; FALSE
# where the old start should run instead, with nothing left behind.
# `dir`: `.ctJuliaServerReady()`'s.
.ctJuliaServerStart <- function(wire, internals, dir) {
  exe <- tryCatch(internals$exe(), error = function(e) NULL)
  if (is.null(exe)) return(FALSE)
  portfile <- tempfile(paste0("juliaPort", Sys.getpid()))
  output <- tempfile("startupOutput")
  on.exit(unlink(c(portfile, output)), add = TRUE)
  opts <- Sys.getenv("JULIACONNECTOR_JULIAOPTS")
  args <- c(shQuote(paste0("--project=", file.path(dir, "env"))),
    if (nzchar(opts)) opts, "--", shQuote(file.path(dir, "serve.jl")), "11980",
    shQuote(portfile), "f")
  started <- tryCatch({
    system2(command = exe, args = args, wait = FALSE, stdout = output, stderr = output,
      env = internals$env())
    TRUE
  }, error = function(e) FALSE)
  if (!started) return(FALSE)
  # JuliaConnectoR's own patience for a server that is starting, and no more
  # for one that has failed: an error in its output ends the wait.
  waited <- 0
  silence <- as.integer(Sys.getenv("JULIACONNECTOR_SILENCE_TIMEOUT", 60))
  quiet <- 0
  size <- 0
  while (file.access(portfile, mode = 4) < 0) {
    Sys.sleep(0.05)
    waited <- waited + 0.05
    now <- suppressWarnings(file.info(output)$size)
    if (!is.na(now) && now > size) {
      size <- now
      quiet <- 0
      text <- tryCatch(readLines(output, warn = FALSE), error = function(e) character())
      if (any(grepl("^ERROR:", text))) return(FALSE)
    } else quiet <- quiet + 0.05
    if (quiet >= silence || waited >= 300) return(FALSE)
  }
  port <- tryCatch(suppressWarnings(as.integer(readLines(portfile, n = 1L, warn = FALSE))),
    error = function(e) NA_integer_)
  if (length(port) != 1L || is.na(port)) return(FALSE)
  # JuliaConnectoR reads the address only as it connects; unset at once, so
  # nothing later takes this session for one someone else started.
  Sys.setenv(JULIACONNECTOR_SERVER = paste0("localhost:", port))
  on.exit(Sys.unsetenv("JULIACONNECTOR_SERVER"), add = TRUE)
  connected <- tryCatch({
    withCallingHandlers(wire$ensure(), message = function(m) {
      if (startsWith(conditionMessage(m), "Connecting to Julia TCP server")) {
        invokeRestart("muffleMessage")
      }
    })
    !is.null(wire$pkgLocal$con)
  }, error = function(e) FALSE)
  if (!connected) {
    # The server serves one client and ends once it has gone; one that never
    # had one would wait for it for good.
    if (!is.null(wire$pkgLocal$con)) {
      suppressWarnings(try(JuliaConnectoR::stopJulia(), silent = TRUE))
    } else {
      try(close(socketConnection("localhost", port, open = "r+b", timeout = 2)),
        silent = TRUE)
    }
  }
  connected
}

# The precompiled server to start this session from: its package's directory
# and JuliaConnectoR's internals, or NULL where there is none to be had, which
# is then not asked again this R session.
.ctJuliaServerPrepare <- function() {
  if (isTRUE(.ct_julia_cache$rserver_off)) return(NULL)
  out <- tryCatch({
    internals <- .ctJuliaServerInternals()
    bin <- .ctJuliaBin()
    dir <- if (!is.null(internals) && !is.null(bin)) .ctJuliaServerReady(internals, bin)
    if (is.null(dir)) NULL else list(internals = internals, dir = dir)
  }, error = function(e) NULL)
  if (is.null(out)) .ct_julia_cache$rserver_off <- TRUE
  out
}
