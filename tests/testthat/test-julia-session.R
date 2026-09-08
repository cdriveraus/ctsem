# Julia session management, and one whole fit from a model and a data frame.
#
# Everything else in the julia test files reaches for a fitted object, so all of
# them share one setup path and none of them exercises it. That is how
# `.ctJuliaSessionRunning()` came to answer "is a session running?" by *starting*
# one -- `JuliaConnectoR::juliaEval()` starts a server rather than reporting on
# it -- with two silent consequences:
#
#   * `ctJuliaSetup(threads = n, force = TRUE)` cleared the session, asked the
#     question, saw the session its own question had just created, and warned
#     that `threads` would have no effect. `force` could not work.
#   * `ctFit(backend = 'julia', cores = n)` sets `JULIA_NUM_THREADS` only when no
#     session is running, so it was never set and every fit ran on one thread
#     whatever `cores` said.
#
# Both are asserted below, along with a plain model-and-data-to-fit check that
# does not reuse a cached fixture.

test_that("asking whether a session is running does not start one", {
  skip_without_julia()
  ctsem:::.ctJuliaClearSession()
  # The regression: this must not be what brings a session into existence.
  expect_false(ctsem:::.ctJuliaSessionRunning())
  expect_false(ctsem:::.ctJuliaSessionRunning())
})

test_that("force restarts the session and threads take effect", {
  skip_without_julia()
  previous <- Sys.getenv("JULIA_NUM_THREADS", unset = NA)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("JULIA_NUM_THREADS")
    else Sys.setenv(JULIA_NUM_THREADS = previous)
  }, add = TRUE)

  # Two threads is enough to prove the setting is applied without asking the
  # test machine for cores it may not have.
  expect_no_warning(suppressMessages(ctJuliaSetup(threads = 2L, force = TRUE)))
  expect_equal(Sys.getenv("JULIA_NUM_THREADS"), "2")
  expect_equal(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")), 2L)
  expect_true(ctsem:::.ctJuliaSessionRunning())

  # Asking for a different count *without* force cannot take effect, and saying
  # so is the point of the warning -- it must fire here and only here.
  expect_warning(suppressMessages(ctJuliaSetup(threads = 3L)), "already running")
  expect_equal(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")), 2L)
})

test_that("a model and a data frame fit end to end", {
  skip_without_julia()
  set.seed(20260827)
  nsubjects <- 12; nobs <- 6
  dat <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- stats::rnorm(1, 1.5, 0.9)
    state <- stats::rnorm(1, 0, 0.5)
    out <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        decay <- exp(-0.4)
        state <- decay * state +
          stats::rnorm(1, 0, sqrt(0.36 / 0.8 * (1 - decay^2)))
      }
      out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = out)
  }))
  # MANIFESTVAR is pinned at the sd the data were generated with rather than
  # left free. Free process noise *and* free measurement error on a single
  # indicator over 12 subjects and 6 occasions does not identify both: mvarY1
  # walks to raw 10.17, where the variance transform is flat to machine
  # precision, and the saturation guard then correctly reports the fit as not
  # converged -- identically on both routes, with pinned or random starting
  # values. That is the guard working, and this file's subject is the session
  # and worker path, not weak identification. The guard has its own coverage in
  # test-julia-binary.R ("a saturated optimum is not reported as converged"),
  # so pinning here loses nothing.
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), T0VAR = matrix(0.5),
    MANIFESTVAR = matrix(0.3), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE

  for (route in c("laplace", "augmented")) {
    fit <- suppressMessages(ctFit(dat, model, backend = "julia",
      intoverpop = if (route == "laplace") "laplace" else TRUE, cores = 1,
      optimcontrol = list(estonly = TRUE)))
    expect_s3_class(fit, "ctJuliaFit")
    expect_true(isTRUE(fit$estimate$converged), label = route)
    expect_false(isTRUE(fit$estimate$stalled), label = route)
    expect_true(is.finite(fit$estimate$loglik), label = route)
    # The population mean of an identity-transformed MANIFESTMEANS is the one
    # parameter whose truth is known here without any transform bookkeeping.
    summ <- suppressWarnings(summary(fit))
    expect_true(is.list(summ) || is.data.frame(summ), label = route)
  }
})

# The engine's subject-chunk ceiling is session-global Julia state. Three call
# sites wrote it and none put it back, so after a `cores=8` fit the session read
# 8 and after a `cores=1` fit it read 1 -- and every later ctKalman(),
# ctExtract() or ctLOO() inherited whichever fit came last. Performance-only,
# but it makes a timing depend on history, which is what makes one impossible to
# reproduce.
test_that("a fit leaves the engine's chunk ceiling as it found it", {
  skip_without_julia()
  ceiling <- function() as.integer(JuliaConnectoR::juliaEval(
    "ContinuousTimeSEM.ctsem_max_chunks().max_chunks"))
  # Put it back on the way out. This test is about session state not leaking,
  # and leaving a ceiling of 3 behind for every later test file would be a poor
  # way to make that point -- it did, and it changed which chunk count the next
  # file's fits ran at.
  original <- ceiling()
  on.exit(JuliaConnectoR::juliaEval(sprintf(
    "ContinuousTimeSEM.ctsem_set_max_chunks!(%d)", original)), add = TRUE)
  JuliaConnectoR::juliaEval("ContinuousTimeSEM.ctsem_set_max_chunks!(3)")
  before <- ceiling()
  expect_equal(before, 3L)

  set.seed(9)
  dat <- do.call(rbind, lapply(1:8, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * .5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    cores = 2, optimcontrol = list(estonly = TRUE))))

  expect_equal(ceiling(), before)
  # And the tuner's answer is still recorded on the fit, which is what the
  # uncertainty phase reads instead of `cores`.
  expect_true(is.numeric(fit$estimate$chunks) || is.integer(fit$estimate$chunks))
  expect_lte(as.integer(fit$estimate$chunks), 2L)

  # A diagnostic must not reconfigure the session either.
  JuliaConnectoR::juliaEval("ContinuousTimeSEM.ctsem_set_max_chunks!(3)")
  suppressWarnings(suppressMessages(try(ctExtract(fit, cores = 2), silent = TRUE)))
  expect_equal(ceiling(), 3L)
})

# `cores` is a ceiling, and `ctsem_tune_chunks!` may measure a much smaller
# chunk count as fastest within it -- rightly, since the subject loop is not
# monotone in the count. That was silent, so `ctFit(cores = 12)` could run on
# two chunks with nothing said. What is tested here is the gate rather than the
# wording: it must not fire when the tuner used the headroom, and it must not
# blame the tuner for a shortfall that was really the session's thread count.
#
# `threads` is passed rather than taken from the session, so the case exercised
# does not depend on what an earlier file left the session at -- the test above
# leaves it at two, where the gate can never fire and this would have skipped.
test_that("a fit says so when the tuner used far fewer chunks than cores allowed", {
  cache <- ctsem:::.ct_julia_cache
  # The 'said once' record is session state like the ceiling above, so it is put
  # back rather than left for the next file to inherit.
  original <- cache$chunks_reported
  on.exit(cache$chunks_reported <- original, add = TRUE)
  said <- function(cores, picked, threads = 12L) {
    seen <- character()
    withCallingHandlers(
      ctsem:::.ctBackendReportChunks(cores, picked, threads = threads),
      message = function(m) {
        seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage")
      })
    any(grepl("^cores = ", seen))
  }
  # Each probe starts from a clean 'said once' record.
  fired <- function(...) { cache$chunks_reported <- NULL; said(...) }

  # Materially short of the ceiling, so worth a line -- and the line names both
  # numbers and where the count it used is recorded.
  cache$chunks_reported <- NULL
  expect_message(
    ctsem:::.ctBackendReportChunks(12L, 2L, threads = 12L),
    "cores = 12.*2 used")

  # The tuner used what it was given, or nearly, or the shortfall is one core:
  # nothing worth interrupting for.
  expect_false(fired(12L, 12L))
  expect_false(fired(12L, 7L))
  expect_false(fired(2L, 1L))
  expect_false(fired(1L, 1L))
  # And nothing to report from a fit whose engine did not say what it used.
  expect_false(fired(12L, NA_integer_))

  # Asking for more cores than the session has threads is a different story --
  # one about `ctJuliaSetup(threads=)`, not about the tuner, which never saw the
  # wider setting and so never rejected it. Reporting it here would misattribute
  # it, and would name a ceiling the tuner had not measured against.
  expect_false(fired(24L, 4L, threads = 4L))
  expect_true(fired(24L, 4L, threads = 24L))

  # Said once, so a simulation study looping a hundred fits gets one line.
  expect_true(fired(12L, 2L))
  expect_false(said(12L, 2L))
  expect_false(said(12L, 2L))
  # A different answer from the tuner is a different thing to say, though.
  expect_true(said(12L, 4L))
})

test_that("the bridge socket is tuned when the session starts", {
  skip_without_julia()
  # What this guards is the finding rather than the speed: JuliaConnectoR sends
  # each message as a run of small writes, which on Linux costs a ~40 ms
  # Nagle/delayed-ACK stall per message unless the socket says otherwise. The
  # tuning is reached through JuliaConnectoR's private communicator object, so
  # the thing most likely to break it is that package renaming something --
  # which would be silent, because every failure path here is a tryCatch.
  suppressMessages(ctJuliaSetup())
  expect_false(is.null(ctsem:::.ctJuliaCommunicator()))

  # 1 when the quickack task is running, 0 when only Nagle was disabled.
  # TCP_QUICKACK is Linux-only, so 0 is the correct answer everywhere else.
  state <- ctsem:::.ctJuliaTuneBridge()
  expect_true(state %in% c(0L, 1L))
  expect_equal(state, if (Sys.info()[["sysname"]] == "Linux") 1L else 0L)
  expect_equal(as.logical(JuliaConnectoR::juliaEval(
    "ContinuousTimeSEM.ctsem_bridge_tuned()")), state == 1L)

  # Idempotent: a second call must not leave a second task re-arming the socket.
  expect_equal(ctsem:::.ctJuliaTuneBridge(), state)

  # And the option turns the whole thing off, for anyone who needs the socket
  # left exactly as JuliaConnectoR opened it.
  withr::with_options(list(ctsem.julia.tunebridge = FALSE), {
    expect_true(is.na(ctsem:::.ctJuliaTuneBridge()))
  })
})

# The warmed sampling worker pool (R/ctBackendWarmWorkers.R) had no test
# coverage at all before this: .ctBackendWarmStop() was defined and never
# called from anywhere in R/ or tests/, so nothing here proved the pool a
# caller starts is ever actually released. `future::nbrOfWorkers()` is what
# the pool itself is built on (a `multisession` plan), so it is what counts
# the workers without needing to reach into `future`'s internals.
test_that("starting and stopping the warmed pool leaves no workers behind", {
  skip_without_julia()
  skip_if_not_installed("future")

  model <- suppressWarnings(ctModel(type = "ct", LAMBDA = diag(1),
    DRIFT = matrix("drift", 1, 1), DIFFUSION = matrix("diffusion", 1, 1),
    MANIFESTVAR = matrix("residual", 1, 1), MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1)))
  dat <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2), Y1 = 0)
  prepared <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))

  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::sequential)
  expect_equal(future::nbrOfWorkers(), 1L)

  handles <- ctsem:::.ctBackendWarmWorkers(prepared, workers = 2L)
  skip_if(is.null(handles), "warming did not start on this machine")
  expect_equal(future::nbrOfWorkers(), 2L)
  expect_false(inherits(future::plan(), "sequential"))

  # The exported release function, not the internal one directly: this is the
  # documented way a caller gets the workers back.
  ctJuliaWorkersStop()
  expect_equal(future::nbrOfWorkers(), 1L)
  expect_true(inherits(future::plan(), "sequential"))
})

test_that("a pool that fails to warm any worker is not left half-started", {
  skip_if_not_installed("future")
  # No live Julia needed: `future::future()` itself is made to fail for every
  # worker, so `future::plan(multisession)` is set (a real, if pointless, pool
  # of background R processes) and then nothing warms in it -- the branch
  # `on.exit(if (!ok) .ctBackendWarmStop(NULL))` exists for. Without it this
  # would return NULL having left the multisession plan and the
  # connections-misuse override in place for a pool that warmed nothing.
  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::sequential)
  previous_option <- getOption("future.connections.onMisuse")

  testthat::local_mocked_bindings(future = function(...) stop("boom"), .package = "future")
  fake <- structure(list(parameter_table = data.frame(parnumber = 1)), class = "ctJuliaModel")
  handles <- ctsem:::.ctBackendWarmWorkers(fake, workers = 2L)

  expect_null(handles)
  expect_true(inherits(future::plan(), "sequential"))
  expect_identical(getOption("future.connections.onMisuse"), previous_option)
})

# `cores` on the julia backend is capped by the Julia session's thread count,
# and Julia fixes that count at process start. So a script that starts the
# engine before its first fit -- which any benchmark harness does, to pay the
# engine load once and outside the timing -- pins every later `cores = n` fit
# to one subject chunk. That ran serially and said nothing, and it invalidated
# two whole benchmark passes: the tell was that the gradient counts for
# `cores = 1` and `cores = 4` came back identical to the digit.
#
# The gate is checked without a session first, since `.ctBackendResolveThreads()`
# takes `threads` for exactly that, and then the whole thing is driven for real:
# a session started at one thread, asked for two, must either get two or say so.
test_that("a fit asking for more cores than the session has threads says so", {
  cache <- ctsem:::.ct_julia_cache
  original <- cache$threads_reported
  on.exit(cache$threads_reported <- original, add = TRUE)
  said <- function(cores, threads, report = TRUE) {
    seen <- character()
    withCallingHandlers(
      ctsem:::.ctBackendResolveThreads(cores, threads = threads,
        report = report),
      message = function(m) {
        seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage")
      })
    seen
  }
  fired <- function(...) {
    cache$threads_reported <- NULL
    any(grepl("^cores = ", said(...)))
  }

  # The case that cost the benchmark, and the line it now prints. It names both
  # numbers, why the shortfall cannot be undone in place, and the exact call
  # that undoes it.
  cache$threads_reported <- NULL
  line <- said(4L, 1L)
  expect_length(line, 1L)
  expect_equal(line, paste0(
    "cores = 4 requested, 1 used: Julia's thread count is fixed at session ",
    "start. ctJuliaSetup(threads = 4, force = TRUE) restarts it with 4.
"))

  # Nothing to say when the session can give what was asked for, when nothing
  # was asked for, or when the count could not be established.
  expect_false(fired(4L, 4L))
  expect_false(fired(4L, 8L))
  expect_false(fired(1L, 1L))
  expect_false(fired(4L, NA_integer_))
  # `report = FALSE` is `fit = FALSE`: a call that only builds a specification
  # has no cores to fall short of, so it takes the thread count it can and says
  # nothing about the one it cannot.
  expect_false(fired(4L, 1L, report = FALSE))
  # Two is the default `cores` and one is Julia's default thread count, so this
  # is the pairing a user meets first. Halving the available parallelism is
  # worth a line as much as quartering it is.
  expect_true(fired(2L, 1L))

  # Said once per pair, so a simulation study looping a hundred fits gets one
  # line -- but a fit that asks for a different number is a different thing to
  # say and is said.
  expect_true(fired(8L, 1L))
  expect_false(any(grepl("^cores = ", said(8L, 1L))))
  expect_true(any(grepl("^cores = ", said(8L, 2L))))
})

test_that("a one-thread session asked for more cores either gets more or says so", {
  skip_without_julia()
  cache <- ctsem:::.ct_julia_cache
  previous_env <- Sys.getenv("JULIA_NUM_THREADS", unset = NA)
  previous_from_cores <- cache$threads_from_cores
  previous_reported <- cache$threads_reported
  nthreads <- function() as.integer(
    JuliaConnectoR::juliaEval("Threads.nthreads()"))
  on.exit({
    cache$threads_from_cores <- previous_from_cores
    cache$threads_reported <- previous_reported
    if (is.na(previous_env)) Sys.unsetenv("JULIA_NUM_THREADS")
    else Sys.setenv(JULIA_NUM_THREADS = previous_env)
  }, add = TRUE)

  set.seed(11)
  dat <- do.call(rbind, lapply(1:8, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * .5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  # Every message, not the first one: `expect_message()` looks at one condition
  # and a fit raises several, so the line under test would be missed by it.
  messages <- function(expr) {
    seen <- character()
    withCallingHandlers(suppressWarnings(force(expr)),
      message = function(m) {
        seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage")
      })
    seen
  }
  onefit <- function() ctFit(dat, model, backend = "julia", cores = 2,
    optimcontrol = list(estonly = TRUE))

  # The harness's own opening move: start the engine, then fit.
  suppressMessages(suppressWarnings(ctJuliaSetup(threads = 1L, force = TRUE)))
  expect_equal(nthreads(), 1L)
  cache$threads_reported <- NULL
  seen <- messages(onefit())
  expect_true(any(grepl("cores = 2 requested, 1 used", seen)),
    label = paste(seen, collapse = " | "))
  # And it really did run on one thread; the message is not decoration.
  expect_equal(nthreads(), 1L)

  # Opted in, the same call gets the threads instead of a line about them.
  cache$threads_reported <- NULL
  seen <- withr::with_options(list(ctsem.julia.restart = TRUE),
    messages(onefit()))
  expect_gte(nthreads(), 2L)
  expect_true(any(grepl("Restarting the Julia session at 2 threads", seen)),
    label = paste(seen, collapse = " | "))
  expect_false(any(grepl("cores = 2 requested", seen)))
})

test_that("a control name the sampler does not read is refused, not ignored", {
  # `control$minEss` on `list(minESS = 100)` is NULL, so the setting was
  # accepted and ignored and the run said nothing -- reported from a real
  # session as a sampler that would not stop at the effective size asked for.
  expect_error(.ctBackendSampleCheckControl(list(minESS = 100)),
    "minESS \\(did you mean minEss\\?\\)")
  expect_error(.ctBackendSampleCheckControl(list(nonsense = 1)),
    "does not read: nonsense")
  # The suggestion is case-insensitive, because case is what goes wrong.
  expect_error(.ctBackendSampleCheckControl(list(MaxDraws = 2)),
    "did you mean maxDraws")
  # An unnamed entry cannot be read either, and says so rather than being
  # silently dropped by the `setdiff`.
  expect_error(.ctBackendSampleCheckControl(list(5)), "must be named")

  # Every name the sampler actually reads passes, and the list is the one the
  # readers use -- a knob added to `.ctBackendSampleControl` and not to
  # `.CT_SAMPLE_CONTROL_NAMES` fails here rather than in a user's script.
  expect_null(.ctBackendSampleCheckControl(list(minEss = 100, warmup = 0L,
    target_accept = 0.9, maxdepth = 12L, max_treedepth = 12L, maxdelta = 900,
    adapt_delta = 0.9, init_scale = 1, adapt_metric = TRUE,
    adapt_effects = FALSE, meanEss = 200, maxDraws = 4000L,
    rhatTarget = 1.01, settleTol = 0, seed = 1L, processes = FALSE,
    callback = function(...) NULL)))
  expect_null(.ctBackendSampleCheckControl(list()))
  expect_null(.ctBackendSampleCheckControl(NULL))
})

test_that("a worker pool left busy by an interrupted run is seen as stale", {
  skip_if_not_installed("future")
  # An interrupted sample leaves its chains running in the pool. R is single
  # threaded, so a busy worker at the moment a fit starts warming can only be
  # an orphan -- and `future::future()` waits for a free worker rather than
  # failing, so warming would queue behind a chain nobody is collecting.
  on.exit(try(future::plan(future::sequential), silent = TRUE), add = TRUE)

  # Says nothing about a sequential plan: the caller replaces that anyway.
  future::plan(future::sequential)
  expect_false(.ctBackendWarmPoolStale(2L))

  future::plan(future::multisession, workers = 2L)
  expect_false(.ctBackendWarmPoolStale(2L))

  orphan <- future::future({ Sys.sleep(30); TRUE }, seed = TRUE)
  # The worker is taken up asynchronously, so this waits for the state rather
  # than assuming it: an assertion here that raced would fail intermittently
  # and be read as the detection not working.
  for (attempt in 1:100) {
    if (future::nbrOfFreeWorkers() < future::nbrOfWorkers()) break
    Sys.sleep(0.1)
  }
  expect_lt(future::nbrOfFreeWorkers(), future::nbrOfWorkers())
  expect_true(.ctBackendWarmPoolStale(2L))

  # And the repair is a repair: releasing the sessions takes the orphan with
  # them, and the pool built next is free.
  .ctBackendWarmStop(NULL)
  expect_s3_class(future::plan("list")[[1L]], "sequential")
  future::plan(future::multisession, workers = 2L)
  expect_false(.ctBackendWarmPoolStale(2L))
  expect_equal(future::value(future::future({ 42L }, seed = TRUE)), 42L)
})
