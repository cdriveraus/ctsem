# Which Julia processes are ctsem's, and stopping the ones that should not run --
#
# A Julia process can outlive the reason it was started. Escape leaves one
# finishing an uncheckpointed step in the background; a crashed R leaves its
# engine running until the engine's next checkpoint notices; an older ctsem
# had no such check at all; sampling workers each hold one of their own. None
# of that is visible from R, and on a shared machine several R sessions' worth
# of them look alike in a task manager: `julia.exe`, nothing else.
#
# So each session ctsem starts is written down: a small file per Julia process
# under R_user_dir("ctsem", "cache"), naming the R process it serves, the host,
# when it started and its thread count, and removed when ctsem ends it or R
# exits. That record is what tells a ctsem engine from any other Julia and an
# orphan from a live session's engine. Processes it cannot account for -- a
# JuliaConnectoR session some other code started, or one from a ctsem before
# the record existed -- are found by their parent: a Julia whose parent is an R
# process is a JuliaConnectoR session, and one whose parent is such a Julia is
# its child (precompiling the engine, typically).
#
# The detail comes from the `ps` package when it is installed. Without it the
# record still says which processes exist and whose they are, and nothing
# more; the command line is no substitute on Windows, where `ps` cannot read
# it and a JuliaConnectoR session shows only the path to julia.exe.

.ctJuliaRegistryDir <- function() {
  getOption("ctsem.julia.registry",
    file.path(tools::R_user_dir("ctsem", which = "cache"), "julia", "processes"))
}

.ctJuliaHost <- function() {
  gsub("[^A-Za-z0-9._-]", "_", Sys.info()[["nodename"]])
}

# Written by .ctJuliaEnsureStarted() once the process id is known. One file per
# process, named by host as well as id, since R_user_dir() can be a network
# home shared by machines whose process ids mean nothing to each other.
.ctJuliaRegister <- function(pid, threads = NA_integer_) {
  if (is.null(pid) || is.na(pid)) return(invisible(FALSE))
  dir <- .ctJuliaRegistryDir()
  ok <- tryCatch({
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    record <- data.frame(julia_pid = as.integer(pid), r_pid = Sys.getpid(),
      host = .ctJuliaHost(), started = format(as.numeric(Sys.time()), nsmall = 3),
      threads = as.integer(threads),
      engine = tryCatch(.ctJuliaEngineVersion(), error = function(e) NA_character_))
    write.dcf(record, .ctJuliaRegistryFile(pid))
    TRUE
  }, error = function(e) FALSE)
  invisible(ok)
}

.ctJuliaRegistryFile <- function(pid) {
  file.path(.ctJuliaRegistryDir(), paste0(.ctJuliaHost(), "-", as.integer(pid), ".dcf"))
}

.ctJuliaUnregister <- function(pid) {
  if (is.null(pid) || !length(pid) || all(is.na(pid))) return(invisible(NULL))
  unlink(.ctJuliaRegistryFile(pid[!is.na(pid)]))
  invisible(NULL)
}

# This host's records, with the ones whose process has gone deleted as read.
.ctJuliaRegistry <- function(info = NULL) {
  dir <- .ctJuliaRegistryDir()
  files <- list.files(dir, pattern = paste0("^", .ctJuliaHost(), "-\\d+\\.dcf$"),
    full.names = TRUE)
  empty <- data.frame(julia_pid = integer(), r_pid = integer(), started = numeric(),
    threads = integer(), engine = character(), stringsAsFactors = FALSE)
  if (!length(files)) return(empty)
  rows <- lapply(files, function(f) tryCatch({
    d <- as.data.frame(read.dcf(f), stringsAsFactors = FALSE)
    data.frame(julia_pid = as.integer(d$julia_pid), r_pid = as.integer(d$r_pid),
      started = as.numeric(d$started), threads = as.integer(d$threads),
      engine = as.character(d$engine), file = f, stringsAsFactors = FALSE)
  }, error = function(e) NULL))
  table <- do.call(rbind, rows)
  if (is.null(table) || !nrow(table)) return(empty)
  if (is.null(info)) info <- .ctProcInfo(table$julia_pid)
  info <- info[match(table$julia_pid, info$pid), , drop = FALSE]
  # A record whose process has gone, or whose id now belongs to a process that
  # is not Julia or that started well before or after the record, is stale:
  # Windows hands process ids out again.
  created <- info$created
  plausible <- !is.na(info$alive) & info$alive &
    grepl("^julia", info$name, ignore.case = TRUE) &
    (is.na(created) | (created <= table$started + 5 & created >= table$started - 3600))
  stale <- !plausible
  if (any(stale)) unlink(table$file[stale])
  table <- table[!stale, setdiff(names(table), "file"), drop = FALSE]
  rownames(table) <- NULL
  table
}

# What can be learned about some processes: alive, name, parent, start time,
# CPU seconds, resident memory, OS threads. With `ps` all of it; without, only
# whether each is running and what its executable is called.
.ctProcInfo <- function(pids, use_ps = requireNamespace("ps", quietly = TRUE)) {
  pids <- unique(as.integer(pids[!is.na(pids)]))
  out <- data.frame(pid = pids, alive = FALSE, name = NA_character_,
    ppid = NA_integer_, parent = NA_character_, created = NA_real_,
    cpu = NA_real_, rss_mb = NA_real_, os_threads = NA_integer_,
    stringsAsFactors = FALSE)
  if (!length(pids)) return(out)
  if (isTRUE(use_ps)) {
    for (i in seq_along(pids)) {
      h <- tryCatch(ps::ps_handle(pids[i]), error = function(e) NULL)
      if (is.null(h) || !isTRUE(tryCatch(ps::ps_is_running(h), error = function(e) FALSE))) next
      safe <- function(expr, default) tryCatch(expr, error = function(e) default)
      out$alive[i] <- TRUE
      out$name[i] <- safe(ps::ps_name(h), NA_character_)
      out$ppid[i] <- safe(as.integer(ps::ps_ppid(h)), NA_integer_)
      out$parent[i] <- safe(ps::ps_name(ps::ps_parent(h)), NA_character_)
      out$created[i] <- safe(as.numeric(ps::ps_create_time(h)), NA_real_)
      out$cpu[i] <- safe(sum(ps::ps_cpu_times(h)[c("user", "system")]), NA_real_)
      out$rss_mb[i] <- safe(ps::ps_memory_info(h)[["rss"]] / 2^20, NA_real_)
      out$os_threads[i] <- safe(as.integer(ps::ps_num_threads(h)), NA_integer_)
    }
    return(out)
  }
  # Every process, not only Julia: these ids are often R sessions, and one
  # missing from a Julia-only list would read as ended -- which is what makes
  # its engine an orphan that ctJuliaKill() stops by default.
  listed <- .ctProcList(pattern = ".", use_ps = FALSE)
  hit <- match(out$pid, listed$pid)
  out$alive <- !is.na(hit)
  out$name <- listed$name[hit]
  out
}

# Running processes matching `pattern`, as (pid, name): by default the Julia
# processes a record does not cover. Through `ps` by id and name only -- 20 ms
# for 500 processes on Windows, where `ps::ps()` takes six seconds collecting
# everything else and `tasklist` most of one; the system's own lister without.
.ctProcList <- function(pattern = "^julia", use_ps = requireNamespace("ps", quietly = TRUE)) {
  empty <- data.frame(pid = integer(), name = character(), stringsAsFactors = FALSE)
  rows <- tryCatch({
    if (isTRUE(use_ps)) {
      pids <- as.integer(ps::ps_pids())
      names <- vapply(pids, function(p) tryCatch(ps::ps_name(ps::ps_handle(p)),
        error = function(e) NA_character_), character(1))
      data.frame(pid = pids, name = names, stringsAsFactors = FALSE)
    } else if (.Platform$OS.type == "windows") {
      lines <- suppressWarnings(system2("tasklist", c("/FO", "CSV", "/NH"),
        stdout = TRUE, stderr = FALSE))
      fields <- lapply(strsplit(lines, "\",\"", fixed = TRUE), function(x) gsub("\"", "", x))
      fields <- Filter(function(x) length(x) >= 2L, fields)
      data.frame(pid = as.integer(vapply(fields, `[`, "", 2L)),
        name = vapply(fields, `[`, "", 1L), stringsAsFactors = FALSE)
    } else {
      lines <- suppressWarnings(system2("ps", c("-A", "-o", "pid=,comm="), stdout = TRUE))
      parts <- regmatches(lines, regexec("^\\s*(\\d+)\\s+(.*)$", lines))
      parts <- Filter(function(x) length(x) == 3L, parts)
      data.frame(pid = as.integer(vapply(parts, `[`, "", 2L)),
        name = basename(vapply(parts, `[`, "", 3L)), stringsAsFactors = FALSE)
    }
  }, error = function(e) empty)
  rows <- rows[!is.na(rows$pid) & grepl(pattern, rows$name, ignore.case = TRUE), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

.ctProcIsR <- function(name) {
  !is.na(name) & grepl("^(R|Rterm|Rscript|Rgui|rsession|R-term)(\\.exe)?$", name,
    ignore.case = TRUE)
}

#' Julia processes started for ctsem
#'
#' Lists the Julia processes on this machine that ctsem started or that
#' belong to them: this R session's engine, the engines of its sampling
#' workers, those of other R sessions, and any left running by an R session
#' that has ended. A Julia process can outlive what it was started for --
#' Escape leaves one finishing an uninterruptible step in the background, and a
#' crashed R session leaves its engine running until the engine notices -- so
#' this is where to look when a machine is busier than it should be.
#' \code{\link{ctJuliaKill}} stops them.
#'
#' ctsem records each Julia session it starts, which is how its processes are
#' told from any other Julia. A JuliaConnectoR session it has no record of --
#' one started by other code, or by a ctsem too old to record it -- is listed
#' too, as is a Julia process started by one of these (precompiling the engine,
#' typically). CPU, memory, thread and parent details need the \pkg{ps}
#' package; without it, only which processes exist and whose they are.
#'
#' @param interval Seconds over which to measure CPU use, which is what the
#'   \code{busy} column reports; \code{0} skips the measurement.
#' @return A data frame, one row per process, with columns \code{pid};
#'   \code{role} (\code{"this session"}, \code{"worker of this session"},
#'   \code{"another R session"}, \code{"orphaned"} for one whose R session has
#'   ended, \code{"JuliaConnectoR, unrecorded"}, or \code{"child of <pid>"});
#'   \code{busy} and \code{cpu_percent} over \code{interval};
#'   \code{cpu_seconds}; \code{memory_mb}; \code{julia_threads}, the Julia
#'   thread count ctsem started it with; \code{os_threads}; \code{r_pid}, the
#'   R process it serves; \code{started}; and \code{note}, which says when
#'   this session's engine is still finishing an interrupted call.
#' @seealso \code{\link{ctJuliaKill}}, \code{\link{ctJuliaStatus}}
#' @examples
#' \donttest{
#' ctJuliaProcesses()
#' }
#' @export
ctJuliaProcesses <- function(interval = 0.5) {
  candidates <- .ctProcList()
  reg <- .ctJuliaRegistry()
  pids <- unique(c(reg$julia_pid, candidates$pid))
  info <- .ctProcInfo(pids)
  info <- info[info$alive, , drop = FALSE]
  reg <- reg[reg$julia_pid %in% info$pid, , drop = FALSE]

  me <- Sys.getpid()
  rparents <- .ctProcInfo(unique(c(reg$r_pid, info$ppid)))
  r_alive <- function(pid) {
    hit <- match(pid, rparents$pid)
    !is.na(hit) & rparents$alive[hit]
  }
  r_parent_of <- function(pid) rparents$ppid[match(pid, rparents$pid)]

  role <- rep(NA_character_, nrow(info))
  r_pid <- rep(NA_integer_, nrow(info))
  threads <- rep(NA_integer_, nrow(info))
  started <- info$created
  for (i in seq_len(nrow(info))) {
    k <- match(info$pid[i], reg$julia_pid)
    if (!is.na(k)) {
      r_pid[i] <- reg$r_pid[k]
      threads[i] <- reg$threads[k]
      if (is.na(started[i])) started[i] <- reg$started[k]
      role[i] <- if (identical(reg$r_pid[k], me)) "this session" else
        if (!r_alive(reg$r_pid[k])) "orphaned" else
        if (isTRUE(r_parent_of(reg$r_pid[k]) == me)) "worker of this session" else
        "another R session"
    } else if (.ctProcIsR(info$parent[i])) {
      r_pid[i] <- info$ppid[i]
      role[i] <- if (identical(info$ppid[i], me)) "this session" else "JuliaConnectoR, unrecorded"
    }
  }
  # A Julia started by one of those: a precompile, most often.
  for (i in which(is.na(role))) {
    parent <- info$ppid[i]
    if (!is.na(parent) && parent %in% info$pid[!is.na(role)]) {
      role[i] <- paste0("child of ", parent)
      r_pid[i] <- r_pid[match(parent, info$pid)]
    }
  }
  keep <- !is.na(role)
  info <- info[keep, , drop = FALSE]
  role <- role[keep]; r_pid <- r_pid[keep]; threads <- threads[keep]; started <- started[keep]

  cpu_percent <- rep(NA_real_, nrow(info))
  if (nrow(info) && isTRUE(interval > 0) && any(!is.na(info$cpu))) {
    t0 <- Sys.time()
    Sys.sleep(interval)
    again <- .ctProcInfo(info$pid)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cpu_percent <- round(100 * (again$cpu[match(info$pid, again$pid)] - info$cpu) / elapsed)
  }
  note <- rep("", nrow(info))
  mine <- role == "this session" & info$pid %in% .ctJuliaOr(.ct_julia_cache$pid, NA_integer_)
  if (!is.null(.ct_julia_cache$inflight)) note[mine] <- "finishing an interrupted call"

  out <- data.frame(pid = info$pid, role = role,
    busy = !is.na(cpu_percent) & cpu_percent >= 10, cpu_percent = cpu_percent,
    cpu_seconds = round(info$cpu, 1), memory_mb = round(info$rss_mb),
    julia_threads = threads, os_threads = info$os_threads, r_pid = r_pid,
    started = as.POSIXct(started, origin = "1970-01-01"), note = note,
    stringsAsFactors = FALSE)
  order <- c("this session", "worker of this session", "orphaned", "another R session",
    "JuliaConnectoR, unrecorded")
  out <- out[order(match(sub("^child of .*", "zz", out$role), c(order, "zz")), out$pid), ,
    drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Stop Julia processes started for ctsem
#'
#' Stops the Julia processes \code{\link{ctJuliaProcesses}} lists. By default
#' only those left by R sessions that have ended, which nothing can be using.
#'
#' @param which Which to stop: \code{"orphaned"} (the default) for processes
#'   whose R session has ended; \code{"session"} for this R session's engine,
#'   which the next julia call replaces with a new one (losing what it had
#'   compiled); \code{"workers"} for this session's sampling workers;
#'   \code{"mine"} for all three; \code{"all"} for every process listed,
#'   including the engines of other running R sessions, which then fail.
#' @param pid Process ids to stop instead of choosing by \code{which}. Only
#'   processes \code{ctJuliaProcesses()} lists can be stopped this way.
#' @param ask Whether to ask before stopping another R session's engine, as
#'   \code{which = "all"} or \code{pid} can.
#' @return The rows of \code{ctJuliaProcesses()} that were stopped,
#'   invisibly.
#' @seealso \code{\link{ctJuliaProcesses}}
#' @examples
#' \donttest{
#' ctJuliaKill()             # engines left by R sessions that have ended
#' ctJuliaKill("session")    # this session's engine, e.g. stuck compiling
#' }
#' @export
ctJuliaKill <- function(which = c("orphaned", "session", "workers", "mine", "all"),
  pid = NULL, ask = interactive()) {
  which <- match.arg(which)
  procs <- ctJuliaProcesses(interval = 0)
  children <- function(sel) {
    parents <- procs$pid[sel]
    sel | procs$role %in% paste0("child of ", parents)
  }
  if (!is.null(pid)) {
    unknown <- setdiff(as.integer(pid), procs$pid)
    if (length(unknown)) {
      message("Not a ctsem Julia process, so left alone: ", paste(unknown, collapse = ", "), ".")
    }
    sel <- procs$pid %in% as.integer(pid)
  } else {
    sel <- switch(which,
      orphaned = procs$role == "orphaned",
      session = procs$role == "this session",
      workers = procs$role == "worker of this session",
      mine = procs$role %in% c("orphaned", "this session", "worker of this session"),
      all = rep(TRUE, nrow(procs)))
  }
  sel <- children(sel)
  targets <- procs[sel, , drop = FALSE]
  if (!nrow(targets)) {
    message("No Julia processes to stop.")
    return(invisible(targets))
  }
  others <- targets$role %in% c("another R session", "JuliaConnectoR, unrecorded")
  if (any(others) && isTRUE(ask) &&
      !isTRUE(utils::askYesNo(paste0("Stop ", sum(others), " Julia process(es) that ",
        "other running R sessions are using? Their work will fail."), default = FALSE))) {
    targets <- targets[!others, , drop = FALSE]
    if (!nrow(targets)) return(invisible(targets))
  }
  # This session's own engine goes through the bridge, so the cache forgets it
  # and the next call starts a new one; its workers through the pool, so
  # nothing tries to reuse them.
  if (any(targets$role == "this session")) .ctJuliaKill()
  if (any(targets$role == "worker of this session")) try(.ctBackendWarmStop(NULL), silent = TRUE)
  stopped <- vapply(targets$pid, .ctProcKill, logical(1))
  .ctJuliaUnregister(targets$pid)
  targets$stopped <- stopped | !vapply(targets$pid, function(p)
    isTRUE(.ctProcInfo(p)$alive), logical(1))
  kinds <- table(sub("^child of .*", "child process", targets$role[targets$stopped]))
  message("Stopped ", sum(targets$stopped), " Julia process",
    if (sum(targets$stopped) == 1L) "" else "es",
    if (length(kinds)) paste0(" (", paste(kinds, names(kinds), collapse = ", "), ")"), ".")
  invisible(targets)
}

# End one process, checking first that the id still names a Julia: `ps` holds
# the process's start time in its handle and will not signal a successor.
.ctProcKill <- function(pid) {
  if (requireNamespace("ps", quietly = TRUE)) {
    h <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
    if (is.null(h)) return(FALSE)
    if (!grepl("^julia", tryCatch(ps::ps_name(h), error = function(e) ""), ignore.case = TRUE)) {
      return(FALSE)
    }
    return(isTRUE(tryCatch({ ps::ps_kill(h); TRUE }, error = function(e) FALSE)))
  }
  info <- .ctProcInfo(pid)
  if (!isTRUE(info$alive) || !grepl("^julia", info$name, ignore.case = TRUE)) return(FALSE)
  isTRUE(tools::pskill(pid, tools::SIGKILL))
}

# Said once per R session, when its first engine starts: engines left running
# by R sessions that have ended. From the record alone, so it costs a few
# file reads and no process listing.
.ctJuliaWarnOrphans <- function() {
  if (isTRUE(.ct_julia_cache$orphans_checked)) return(invisible(NULL))
  .ct_julia_cache$orphans_checked <- TRUE
  reg <- tryCatch(.ctJuliaRegistry(), error = function(e) NULL)
  if (is.null(reg) || !nrow(reg)) return(invisible(NULL))
  reg <- reg[reg$r_pid != Sys.getpid(), , drop = FALSE]
  if (!nrow(reg)) return(invisible(NULL))
  parents <- .ctProcInfo(reg$r_pid)
  gone <- !parents$alive[match(reg$r_pid, parents$pid)]
  n <- sum(gone, na.rm = TRUE)
  if (n) {
    message("  ", n, " Julia process", if (n == 1L) " is" else "es are",
      " still running for R sessions that have ended; ctJuliaProcesses() lists ",
      "them, ctJuliaKill() stops them.")
  }
  invisible(NULL)
}
