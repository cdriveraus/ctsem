# ctSample(): Hamiltonian sampling of a julia Laplace fit.
#
# The statistics of the sampler are tested in the engine suite, against
# posteriors known in closed form. What is tested here is the R side of it: that
# a sampled fit is a `ctJuliaFit` the existing summary machinery can read, that
# the diagnostics arrive intact, and that the failure modes are refused rather
# than half-performed.
#
# Deliberately small and short. A sample that mixes well enough to assert
# R-hat on takes minutes, which does not belong in a unit suite; the engine
# tests carry the correctness argument.

.sample_fixture <- function(nsub = 12L, tp = 6L, seed = 31L) {
  set.seed(seed)
  data <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    intercept <- stats::rnorm(1, 0, 0.8)
    state <- stats::rnorm(1, 0, 0.5)
    y <- numeric(tp)
    for (t in seq_len(tp)) {
      state <- 0.75 * state + stats::rnorm(1, 0, 0.4)
      y[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(tp) - 1, Y1 = y)
  }))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[match(TRUE, model$pars$matrix == "MANIFESTMEANS")] <- TRUE
  suppressWarnings(suppressMessages(ctFit(data, model, backend = "julia",
    cores = 1, intoverpop = "laplace", priors = TRUE,
    optimcontrol = list(finishsamples = 20))))
}

test_that("a sampled fit carries draws the summary machinery can read", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  sampled <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 80, draws = 80, cores = 1)))

  expect_s3_class(sampled, "ctJuliaFit")
  # And as a fit, not as a model spec: "ctFitModel" marks an unfitted handle,
  # and the optimised path gives its fits "ctFit".
  expect_s3_class(sampled, "ctFit")
  # The draws are where an optimised fit's normal-approximation draws live, so
  # everything downstream reads them without knowing which produced them.
  expect_equal(dim(sampled$estimate$rawposterior), c(160L, npar))
  expect_true(all(is.finite(sampled$estimate$rawposterior)))
  expect_equal(colnames(sampled$estimate$rawposterior),
    ctsem:::.ctBackendRawParameterNames(fit, npar))

  # The point estimate becomes the posterior mean, and the Laplace estimate it
  # started from is kept rather than overwritten.
  expect_equal(sampled$estimate$raw,
    as.numeric(colMeans(sampled$estimate$rawposterior)))
  expect_equal(sampled$estimate$laplace_raw, fit$estimate$raw)
  expect_equal(sampled$estimate$se,
    sqrt(diag(stats::cov(sampled$estimate$rawposterior))))

  # And the cached constrained draws describe the new draws, not the old ones.
  expect_false(is.null(sampled$transformedpars))
  expect_equal(dim(sampled$transformedpars$samples),
    dim(sampled$estimate$rawposterior))
  expect_no_error(suppressWarnings(summary(sampled)))
})

test_that("the diagnostics come back per parameter and per chain", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  sampled <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 80, draws = 80, cores = 1)))
  diagnostics <- sampled$sample

  expect_s3_class(diagnostics, "ctSampleDiagnostics")
  expect_length(diagnostics$rhat, npar)
  expect_length(diagnostics$ess, npar)
  expect_equal(names(diagnostics$rhat), colnames(sampled$estimate$rawposterior))
  expect_true(all(diagnostics$rhat > 0.99, na.rm = TRUE))
  expect_true(all(diagnostics$ess > 0, na.rm = TRUE))
  expect_length(diagnostics$stepsize, 2L)
  expect_true(all(diagnostics$stepsize > 0))
  expect_length(diagnostics$ebfmi, 2L)
  expect_true(diagnostics$divergent >= 0L)
  expect_length(diagnostics$accept, 160L)
  # The Laplace estimate the chains started from.
  expect_equal(diagnostics$start, fit$estimate$raw)
  expect_output(print(diagnostics), "ctsem Hamiltonian sample")

  # Everything above is asserted at chains > 1 on purpose, because that is the
  # branch that runs each chain in its own process, and it is the branch whose
  # result-assembly tail was once a separate copy that had lost the class, the
  # start, the constrained draws and the effect summaries. The summaries were
  # only ever checked at chains = 1, which is the one count that never reaches
  # the process path -- so a fit could return no random effects at the default
  # setting and the suite stayed green.
  #
  # `processes` records which path ran rather than gating the assertions on it:
  # these hold either way, and a test that skipped when the path was not taken
  # would go quiet in exactly the case worth watching. Note that under
  # `load_all()` a worker cannot load the development tree, so the path falls
  # back to this session and this runs against the in-process tail; an installed
  # library, as `R CMD check` uses, exercises the other one.
  expect_true(is.logical(diagnostics$processes) &&
    length(diagnostics$processes) == 1L)
  neffects <- length(fit$model_spec$subject_starts) *
    length(fit$model_spec$laplace$re_index)
  expect_length(diagnostics$effect_mean, neffects)
  expect_length(diagnostics$effect_sd, neffects)
  expect_true(all(diagnostics$effect_sd > 0))
})

test_that("the effects come back summarised, or in full when asked for", {
  skip_without_julia()
  fit <- .sample_fixture()
  neffects <- length(fit$model_spec$subject_starts) *
    length(fit$model_spec$laplace$re_index)

  # By default only the summary crosses the bridge: the draws themselves are
  # nsubjects x neffects x chains x draws numbers, and the bridge moves about
  # 1 MB/s.
  lean <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 1, warmup = 60, draws = 60, cores = 1)))
  expect_null(lean$sample$effects)
  expect_length(lean$sample$effect_mean, neffects)
  expect_length(lean$sample$effect_sd, neffects)
  expect_true(all(lean$sample$effect_sd > 0))

  full <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 1, warmup = 60, draws = 60, cores = 1,
      saveEffects = TRUE)))
  expect_equal(dim(full$sample$effects), c(60L, neffects))
  expect_true(all(is.finite(full$sample$effects)))
})

test_that("a sampled fit keeps the exact Hessian it was built from", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)

  # The sampler asks the engine for the exact Hessian at the Laplace estimate,
  # and it is not only the metric's starting point: it is what the fit carries
  # as `uncertainty$hessian` and what the identifiability report reads. The
  # call site used to hand `.ctBackendHessian()` an unclassed list, which the
  # objective lookup rejects, so the request errored before the engine saw it
  # -- and because that failure is caught and warned about, every sampled fit
  # said the engine could not differentiate its gradient and then carried no
  # Hessian at all. The warning is asserted against because it is the only
  # thing that was ever visible.
  # Collected rather than `expect_no_warning`d, because a 40-draw chain warns
  # about its own R-hat and effective sample size and those are expected here.
  warnings <- character()
  sampled <- withCallingHandlers(suppressMessages(
    ctSample(fit, chains = 1, warmup = 40, draws = 40, cores = 1)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("could not differentiate", warnings)))

  hessian <- sampled$uncertainty$hessian
  expect_true(is.matrix(hessian))
  expect_equal(dim(hessian), c(npar, npar))
  expect_true(all(is.finite(hessian)))
  # Exact rather than differenced, so the two triangles agree exactly.
  expect_identical(hessian, t(hessian))

  # And the identifiability report is the Hessian's, not the empty one that a
  # NULL silently produces.
  expect_false(is.null(sampled$identifiability))
  expect_true(is.finite(sampled$identifiability$condition))
})

test_that("processes = TRUE reproduces the in-process draws to numerical noise", {
  skip_without_julia()
  fit <- .sample_fixture()

  # `.ctBackendSampleProcesses`'s own comment: `seed + k - 1` makes each
  # worker draw the stream the in-process chain would have used, so with
  # `warmup = 0` -- the first kept draw sitting as near the shared start as
  # the sampler ever gets -- the two routes should agree to numerical noise,
  # not bit for bit: measured there at 6.2e-10 on the first draw. A normal
  # warmup lets NUTS's chaos carry the difference to order 1 within a few
  # dozen transitions, which is why this stays at `warmup = 0`.
  inprocess <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 0, draws = 3, cores = 2, seed = 777,
      processes = FALSE)))
  viaprocess <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 0, draws = 3, cores = 2, seed = 777,
      processes = TRUE)))

  expect_false(isTRUE(inprocess$sample$processes))
  # If the workers could not be used -- in particular, a `future` worker
  # cannot load an uninstalled development tree, so this is the expected
  # outcome under `devtools::load_all()` -- sampling silently falls back to
  # the in-process path, and comparing that fallback against itself would
  # trivially "pass" without checking anything. Skip rather than pass
  # silently; an installed library, as `R CMD check` uses, takes the branch
  # this test is for.
  skip_if_not(isTRUE(viaprocess$sample$processes),
    "processes = TRUE fell back to in-process sampling in this session")

  expect_equal(dim(viaprocess$estimate$rawposterior),
    dim(inprocess$estimate$rawposterior))
  expect_equal(viaprocess$estimate$rawposterior,
    inprocess$estimate$rawposterior, tolerance = 1e-6)
})

test_that("ctSample refuses what it cannot sample", {
  skip_without_julia()
  expect_error(ctSample(list()), "ctFit\\(backend='julia'\\)")

  # The augmented route carries the effects in the state, so there is no
  # separate posterior over them and no Laplace curvature to build a metric
  # from. Saying so beats sampling something else.
  set.seed(4)
  data <- do.call(rbind, lapply(1:8, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * 0.5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  augmented <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, optimcontrol = list(estonly = TRUE))))
  expect_error(ctSample(augmented), "intoverpop='laplace'")
})

test_that("ctOptimUncertainty() refuses a sampled julia fit instead of silently discarding its posterior", {
  skip_without_julia()
  # ctOptimUncertainty()'s whole premise is a point estimate plus curvature.
  # Before this guard, handing it a sampled fit treated fit$estimate$raw (the
  # posterior mean here, not a mode) as that point estimate and overwrote
  # fit$estimate$rawposterior -- the real draws -- with fresh Gaussian draws
  # built from the curvature there. Wrong in a way nothing downstream would
  # notice, since the replacement is the same shape and a plausible size; see
  # review/J7-sampled-fit-support.md.
  fit <- .sample_fixture()
  sampled <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 40, draws = 40, cores = 1)))
  expect_false(is.null(sampled$sample))
  before <- sampled$estimate$rawposterior

  err <- tryCatch({
    ctOptimUncertainty(sampled, uncertainty = "hessian", finishsamples = 20, cores = 1)
    NA_character_
  }, error = function(e) conditionMessage(e))
  expect_false(is.na(err))
  expect_match(err, "sampled", fixed = TRUE)
  expect_match(err, "rawposterior", fixed = TRUE)

  # The refusal happens before anything is touched.
  expect_identical(sampled$estimate$rawposterior, before)

  # The Laplace fit ctSample() started from is a genuinely optimized fit
  # (fit$sample is NULL there) -- unaffected by the new check, still works.
  expect_false(is.null(fit$estimate$raw))
  out <- ctOptimUncertainty(fit, uncertainty = "hessian", finishsamples = 20, cores = 1)
  expect_s3_class(out, "ctJuliaFit")
})

# ---------------------------------------------------------------------------
# Convergence diagnostics, reported the way the stan backend reports them.
#
# The case behind these: a julia sampling run on intensive-longitudinal data
# returned in five minutes with two chains 384,000 log units apart and a
# parameter whose true value is -0.5 estimated at 29.1. Nothing a reader would
# meet said so -- `summary()` had `mean`, `sd` and quantiles where a `ctStanFit`
# has `n_eff` and `Rhat` too, and the fit-time warnings are gone the moment a
# script wraps its call in `suppressWarnings()`, which every batch script does.
#
# The detection itself is tested on constructed draws rather than on a fit,
# because a test that has to make a real sampler fail on cue is a flaky test.
# The fit tests below check that the columns arrive, on the right scale, and
# only for a genuine sample.

test_that("split R-hat and effective size separate agreeing chains from disagreeing ones", {
  set.seed(11)
  # 2 chains x 400 draws, stacked chain-major, which is the layout
  # `estimate$rawposterior` carries.
  healthy <- cbind(a = stats::rnorm(800), b = stats::rnorm(800))
  good <- ctsem:::.ctBackendDrawDiagnostics(healthy, chains = 2L)
  expect_equal(rownames(good), c("a", "b"))
  expect_true(all(good$Rhat < 1.01))
  expect_true(all(good$n_eff > 100))

  # The failure this exists for: chain 2 somewhere else entirely.
  broken <- healthy
  broken[401:800, "a"] <- broken[401:800, "a"] + 5000
  bad <- ctsem:::.ctBackendDrawDiagnostics(broken, chains = 2L)
  expect_gt(bad["a", "Rhat"], 1.01)
  expect_lt(bad["a", "n_eff"], 10)
  # And it is the parameter that moved, not both.
  expect_lt(bad["b", "Rhat"], 1.01)

  # Nothing to diagnose is NULL rather than a number: an optimised fit's draws
  # come from a covariance fitted at the mode and have no chains at all.
  expect_null(ctsem:::.ctBackendDrawDiagnostics(healthy, chains = NULL))
  expect_null(ctsem:::.ctBackendDrawDiagnostics(healthy[1:3, , drop = FALSE], chains = 1L))
})

test_that("the summary's opening line says whether the chains agreed", {
  set.seed(12)
  healthy <- cbind(a = stats::rnorm(800), b = stats::rnorm(800))
  broken <- healthy
  broken[401:800, "a"] <- broken[401:800, "a"] + 5000

  note <- function(values, divergent = 0L) {
    fit <- list(estimate = list(rawposterior = values),
      sample = list(chains = 2L, draws = 400L, divergent = divergent))
    sections <- list(popmeans = ctsem:::.ctBackendSampleSummary(values, chains = 2L))
    ctsem:::.ctBackendSampleNote(fit, sections, 2L)
  }

  good <- note(healthy)
  expect_match(good, "2 chains x 400 draws")
  expect_match(good, "Worst R-hat")
  expect_false(grepl("have not converged", good, fixed = TRUE))

  bad <- note(broken)
  expect_match(bad, "have not converged", fixed = TRUE)
  expect_match(bad, "(a)", fixed = TRUE)

  # A divergence alone is enough, whatever R-hat says.
  expect_match(note(healthy, divergent = 3L), "3 divergent transitions", fixed = TRUE)
  expect_match(note(healthy, divergent = 3L), "have not converged", fixed = TRUE)
})

test_that("a sampled fit reports n_eff and Rhat where a ctStanFit does, and an optimised one does not", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  sampled <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 120, draws = 120, cores = 1)))
  total <- 2L * 120L

  summarised <- suppressWarnings(summary(sampled))
  expect_true(all(c("n_eff", "Rhat") %in% names(summarised$popmeans)))
  # Stan's column order: the moments, the interval, then the diagnostics.
  expect_equal(names(summarised$popmeans),
    c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat"))
  # Not `npar` rows: the individually varying parameter is reported under
  # popsd rather than popmeans, exactly as on the stan backend.
  expect_gt(nrow(summarised$popmeans), 0L)
  expect_lt(nrow(summarised$popmeans), npar + 1L)
  expect_true(all(is.finite(summarised$popmeans$Rhat)))
  # Sane rather than converged: a 120-draw run is not asked to mix, only to
  # produce numbers that are numbers.
  expect_true(all(summarised$popmeans$Rhat > 0.9))
  expect_true(all(summarised$popmeans$n_eff > 0))
  expect_true(all(summarised$popmeans$n_eff <= total * 1.5))
  expect_true(all(c("n_eff", "Rhat") %in% names(summarised$popsd)))

  # The verdict is on the object, not only in a warning that a batch script
  # suppresses, and it is a decision rather than a table to read.
  expect_true(is.logical(sampled$sample$converged))
  expect_length(sampled$sample$converged, 1L)
  expect_true(is.character(sampled$sample$diagnosis))
  expect_output(print(sampled), "chains converged")

  # And the summary opens with it.
  expect_true(is.character(summarised$sampleNote))
  expect_match(summarised$sampleNote, "2 chains x 120 draws", fixed = TRUE)
  expect_identical(names(summarised)[1L], "sampleNote")
  # A sampled fit did not run an uncertainty pass, and used to say it had.
  expect_match(summarised$uncertaintyNote, "Hamiltonian", fixed = TRUE)
  expect_false(grepl("ctOptimUncertainty", summarised$uncertaintyNote, fixed = TRUE))

  # The optimised fit it started from has draws too -- from a covariance fitted
  # at the mode -- and there is no between-chain variance for those. Same
  # treatment as summary.ctStanFit gives an optimised stan fit.
  optimised <- suppressWarnings(summary(fit))
  expect_false(any(c("n_eff", "Rhat") %in% names(optimised$popmeans)))
  expect_null(optimised$sampleNote)
})

test_that("a run too short to mix says so rather than returning quietly", {
  skip_without_julia()
  fit <- .sample_fixture()
  # Twelve draws off a twelve-subject Laplace fit: too few to adapt and far too
  # few to mix. The assertion is on the verdict, not on a threshold being
  # crossed by a particular margin.
  broken <- suppressWarnings(suppressMessages(
    ctSample(fit, chains = 2, warmup = 12, draws = 12, cores = 1)))
  expect_false(isTRUE(broken$sample$converged))
  expect_gt(length(broken$sample$diagnosis), 0L)
  note <- suppressWarnings(summary(broken))$sampleNote
  expect_match(note, "have not converged|Too few effective draws")
})

test_that("the process-path progress line holds every chain on one line", {
  info <- function(phase, iteration, total, logp, divergent = 0L) {
    list(phase = phase, iteration = iteration, total = total, logp = logp,
      divergent = divergent)
  }
  two <- list(info("warmup", 8L, 500L, -2321546.15),
    info("warmup", 4L, 500L, -29733.03))
  line <- .ctBackendProcessLine(two, 2L, c(200, 200), width = 100L)

  # One line, whatever the chain count: it is overwritten in place, and a
  # carriage return only returns to the start of the last visual row.
  expect_length(line, 1L)
  expect_false(grepl("\n", line, fixed = TRUE))
  expect_lte(nchar(line), 99L)
  # The phase and the target are said once for chains that agree on them; the
  # iteration and the log posterior are per chain, because a single chain stuck
  # in a bad region is the failure this line exists to reveal.
  expect_match(line, "warmup 8, 4/500", fixed = TRUE)
  expect_match(line, "-2.32155e+06", fixed = TRUE)
  expect_match(line, "-29733", fixed = TRUE)
  # And the estimate is the slowest chain's: chain 2 has 496 iterations left at
  # 4 per 200s, so it decides when the run ends.
  expect_match(line, "6h 53m at this rate", fixed = TRUE)

  # Divergences appear once there are any, and not before.
  expect_false(grepl("div", line, fixed = TRUE))
  expect_match(.ctBackendProcessLine(list(info("warmup", 8L, 500L, -12.5, 3L),
    info("warmup", 9L, 500L, -12.1, 0L)), 2L, c(10, 10), width = 100L),
    "div 3, 0", fixed = TRUE)

  # Chains in different phases keep their own group rather than being merged
  # into one count that means neither.
  mixed <- .ctBackendProcessLine(list(info("sampling", 12L, 500L, -12.5),
    info("warmup", 498L, 500L, -12.1)), 2L, c(10, 10), width = 100L)
  expect_match(mixed, "sampling 12/500", fixed = TRUE)
  expect_match(mixed, "warmup 498/500", fixed = TRUE)

  # A chain that has not written yet is absent, not an NA, and no chain having
  # written is nothing to print rather than an empty line.
  expect_match(.ctBackendProcessLine(list(info("warmup", 5L, 500L, -12.5), NULL),
    2L, c(10, 10), width = 100L), "warmup 5/500", fixed = TRUE)
  expect_null(.ctBackendProcessLine(list(NULL, NULL), 2L, c(10, 10)))

  # The narrow case is the one that matters, because 80 columns is the default
  # and four chains do not fit in it with everything on. What has to survive is
  # the log posteriors -- truncation used to take exactly those.
  four <- lapply(1:4, function(k) info("warmup", 300L + k, 500L, -5300 + k, 0L))
  narrow <- .ctBackendProcessLine(four, 4L, rep(600, 4), width = 80L)
  expect_lte(nchar(narrow), 79L)
  expect_match(narrow, "logp -5299, -5298, -5297, -5296", fixed = TRUE)
})

test_that("an effective-size target turns the draw count into a budget", {
  skip_without_julia()
  fit <- .sample_fixture()

  # `minESS` used to do nothing at all without `maxDraws`: the engine extends
  # towards a budget that defaulted to exactly the draws asked for, so there
  # was nothing to extend into and the target could only be reported after the
  # fact. The count asked for is the budget now, so a target met early stops
  # the run -- and can only ever shorten it.
  #
  # `rhatTarget` is raised out of the way because this fixture is twelve
  # subjects and will not reach 1.01; what is under test is the stopping rule,
  # not whether a small fixture mixes.
  met <- suppressWarnings(suppressMessages(ctSample(fit, chains = 2,
    warmup = 100, draws = 200, cores = 1, processes = FALSE,
    sampleControl = list(minESS = 0.5, rhatTarget = 100))))
  expect_lt(met$sample$draws, 200L)
  expect_gt(met$sample$draws, 0L)

  # `minESS = 0` is the off switch, and then it takes exactly what it was
  # asked for. Needed explicitly now that 200 is the default target -- a
  # default fit stops early, which is the point of it.
  full <- suppressWarnings(suppressMessages(ctSample(fit, chains = 2,
    warmup = 100, draws = 200, cores = 1, processes = FALSE,
    sampleControl = list(minESS = 0))))
  expect_equal(full$sample$draws, 200L)

  # And the default is a target rather than an instruction: the twelve-subject
  # fixture will not reach min ESS 200 in 200 draws, so the budget is spent in
  # full and nothing is lost by the default being on.
  default <- suppressWarnings(suppressMessages(ctSample(fit, chains = 2,
    warmup = 100, draws = 200, cores = 1, processes = FALSE)))
  expect_lte(default$sample$draws, 200L)

  # The old spelling is refused with the new one named, rather than dropped by
  # `$` and silently ignored.
  expect_error(ctSample(fit, sampleControl = list(minEss = 100)),
    "did you mean minESS")
})

test_that("a fixed stepsize is what every chain starts from", {
  skip_without_julia()
  fit <- .sample_fixture()

  # With no warmup there is nothing to move the step size from where it began,
  # so this is the one setting where the initial value *is* the value -- and
  # where each chain estimating its own showed up as chains that behaved
  # differently for no reason.
  fixed <- suppressWarnings(suppressMessages(ctSample(fit, chains = 3,
    warmup = 0, draws = 20, cores = 1, processes = FALSE,
    sampleControl = list(stepsize = 0.05))))
  expect_equal(fixed$sample$stepsize, rep(0.05, 3L))

  # Left unset, the chains estimate their own and need not agree; what is
  # asserted is only that the setting is not silently ignored.
  free <- suppressWarnings(suppressMessages(ctSample(fit, chains = 3,
    warmup = 0, draws = 20, cores = 1, processes = FALSE)))
  expect_false(isTRUE(all.equal(free$sample$stepsize, rep(0.05, 3L))))

  # And with warmup it is a starting point rather than the answer.
  moved <- suppressWarnings(suppressMessages(ctSample(fit, chains = 2,
    warmup = 60, draws = 20, cores = 1, processes = FALSE,
    sampleControl = list(stepsize = 0.05))))
  expect_false(isTRUE(all.equal(moved$sample$stepsize, rep(0.05, 2L))))
})

# --- the posterior predictive uses the effect draws, not conditional modes ----
#
# `ctsem_generate(laplace, ...)` puts each subject at its conditional mode,
# solved from that subject's own data. For a fit that integrated the random
# effects out that is the quantity to use. For a fit that *sampled* them it is
# not: the draws exist, and generating at modes conditions the predictive on a
# point estimate of each effect instead of integrating over its posterior.
#
# Asserted on the mechanism rather than on the generated data's statistics,
# because the two routes are not far apart in distribution -- with
# fullposterior=TRUE the modes are re-solved at each parameter draw, so they
# move with it, and the aligned correlation between a subject's generated mean
# and its own effect draw measured 0.378 through the draws against 0.241
# through the modes. A difference in the right direction, and far too small a
# gap to test on. The per-subject exactness below is decisive instead.
test_that("generation at given effects is exact and per-subject", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctSample(.sample_fixture(),
    chains = 1, warmup = 60, draws = 60, cores = 1, saveEffects = TRUE)))
  expect_false(is.null(fit$sample$effects))

  raw <- fit$estimate$raw
  nrows <- length(fit$model_spec$times)
  nsub <- length(fit$model_spec$subject_starts)
  rowid <- rep(seq_len(nsub), each = nrows / nsub)
  set.seed(1)
  base <- matrix(stats::rnorm(nrows), 1, nrows)
  effects <- as.numeric(fit$sample$effects[1, ])

  first <- .ctBackendGenerate(fit, raw, base, effects = effects)
  again <- .ctBackendGenerate(fit, raw, base, effects = effects)
  # The effects are data, not a draw: the same ones reproduce the dataset.
  expect_identical(as.numeric(first$Y), as.numeric(again$Y))

  # One subject's effect moves that subject and nothing else. The filter runs
  # per subject at that subject's own parameter vector, so the others are not
  # merely close -- they are the same numbers. Anything that reached them would
  # mean an effect had been applied to the wrong subject, which is the failure
  # this is here to catch.
  moved <- effects
  moved[1] <- moved[1] + 2
  shifted <- .ctBackendGenerate(fit, raw, base, effects = moved)
  y0 <- as.numeric(first$Y)
  y1 <- as.numeric(shifted$Y)
  expect_gt(abs(mean(y1[rowid == 1]) - mean(y0[rowid == 1])), 0.1)
  expect_identical(y1[rowid != 1], y0[rowid != 1])

  # And the effects are actually being used: modes are a different answer.
  modes <- .ctBackendGenerate(fit, raw, base)
  expect_false(isTRUE(all.equal(as.numeric(modes$Y), y0)))
})

test_that("a sampled fit without saved effects says it fell back to modes", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctSample(.sample_fixture(),
    chains = 1, warmup = 60, draws = 60, cores = 1, saveEffects = TRUE)))
  # Silence is the failure mode: the draws are summarised by default, so
  # without a message a user asking for a posterior predictive would get one
  # conditioned on point estimates and no way to notice.
  stripped <- fit
  stripped$sample$effects <- NULL
  expect_message(ctGenerateFromFit(stripped, nsamples = 3,
    fullposterior = TRUE, cores = 1), "saveEffects=TRUE")
  # With them, nothing to report.
  expect_no_message(ctGenerateFromFit(fit, nsamples = 3,
    fullposterior = TRUE, cores = 1))
})

# Which route generates the states, as a choice rather than a consequence.
#
# What a fit contributes to generated data is its sampled individual
# differences. The states are always regenerated -- nothing a fit sampled is
# carried over -- and `intoverstates` says how: through the filter's
# one-step-ahead predictive, or by drawing the trajectory from the process at
# each subject's own parameters and then each observation given its state.
# `'fit'` mirrors what the fit did, which is the default because it keeps the
# identity that a generated dataset's llrow is the likelihood reported while
# generating it.
test_that("ctGenerateFromFit can resample the trajectory or mirror the fit", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctSample(.sample_fixture(),
    chains = 1, warmup = 60, draws = 60, cores = 1, saveEffects = TRUE)))
  expect_true(isTRUE(fit$args$resolved$intoverstates))

  generate <- function(io) {
    set.seed(11)
    suppressMessages(ctGenerateFromFit(fit, nsamples = 6, fullposterior = TRUE,
      cores = 1, intoverstates = io))$generated$Y
  }
  mirror <- generate("fit")
  forced <- generate(TRUE)
  resampled <- generate(FALSE)

  # This fit integrated its states, so mirroring it is the filter route.
  expect_equal(mirror, forced)
  expect_false(isTRUE(all.equal(mirror, resampled)))
  expect_true(all(is.finite(resampled)))
  # Not merely finite: a route that had lost the process or the measurement
  # model would still be finite and would not be on this scale.
  expect_equal(stats::sd(resampled), stats::sd(mirror), tolerance = 0.25)

  # The state-explicit generator and the innovation count both used to throw
  # for a random-effect objective, which is what made this route unreachable
  # from a fit that had any.
  expect_gt(.ctBackendStateDimension(fit), 0L)
})

test_that("the resampled trajectory is drawn at each subject's own parameters", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctSample(.sample_fixture(),
    chains = 1, warmup = 60, draws = 60, cores = 1, saveEffects = TRUE)))
  raw <- fit$estimate$raw
  nrows <- length(fit$model_spec$times)
  nsub <- length(fit$model_spec$subject_starts)
  rowid <- rep(seq_len(nsub), each = nrows / nsub)
  set.seed(2)
  base <- matrix(stats::rnorm(nrows), 1, nrows)
  set.seed(3)
  z <- stats::rnorm(.ctBackendStateDimension(fit))
  effects <- as.numeric(fit$sample$effects[1, ])
  moved <- effects
  moved[1] <- moved[1] + 2

  y0 <- as.numeric(.ctBackendGenerateStates(fit, raw, z, base, effects = effects)$Y)
  y1 <- as.numeric(.ctBackendGenerateStates(fit, raw, z, base, effects = moved)$Y)
  expect_gt(abs(mean(y1[rowid == 1]) - mean(y0[rowid == 1])), 0.1)
  # Bit-identical, as on the filter route: an effect reaching another subject
  # would mean it had been applied to the wrong one.
  expect_identical(y1[rowid != 1], y0[rowid != 1])
})

test_that("intoverstates is refused on the stan path rather than ignored", {
  # An argument that means something on one backend and nothing on the other is
  # the shape of mistake ctsem keeps paying for, so it is refused by name.
  expect_error(ctGenerateFromFit(ctstantestfit, nsamples = 2,
    intoverstates = FALSE), "backend='julia'", fixed = TRUE)
  expect_error(ctGenerateFromFit(ctstantestfit, nsamples = 2,
    intoverstates = TRUE), "backend='julia'", fixed = TRUE)
})
