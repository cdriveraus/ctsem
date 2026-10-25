# ctFitUncertainty(fit, uncertainty = 'sample'): sampling a julia backend fit,
# by SAEM's kernel on the joint posterior and NUTS on the marginal.
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
#
# `ctFitUncertainty(fit, 'sample')` on the `.sample_fixture()` below -- an
# maximum-likelihood fit that named `intoverpop='laplace'` -- defaults to
# the Laplace marginal; had 'auto' chosen the route it would be the joint
# posterior over parameters and random effects (`.ctBackendSampleRoute()`).
# Tests that are about the random effects or SAEM ask for
# `intoverpop = 'none'` explicitly, so they say what they need whatever the
# default; everything else here is generic sampler mechanics that holds under
# either posterior.
#
# The separate exported sampling function this file used to call was removed
# (it was julia-only and never released); every call below that used to reach
# it now reaches `ctFitUncertainty()` instead, with the same settings folded
# into `control` -- `chains`, `warmup`, `draws`, `seed`, `saveEffects`,
# `processes` and `intoverpop` are all entries of one list now, where that
# function had them split across its own arguments and a separate list.

# One fit, shared by every test that asks for it: `laplace_fixture()`
# (helper-julia.R), which test-backend-generate.R generates from too. Each test
# here samples from it or asks it something, and none writes into it. It was
# refitted sixteen times.
.sample_fixture <- function() laplace_fixture()

test_that("a sampled fit carries draws the summary machinery can read", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  # intoverpop='none': the placement asserted below is SAEM's, which SAEM's
  # kernel takes on the joint target; the fixture names 'laplace', whose own
  # target is marginal. Named, because on this fixture, whose effect the filter
  # integrates, the joint default is the map sampler, placed at the fit's
  # estimate (tested below).
  sampled <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 2, warmup = 80, draws = 80, intoverpop = "none",
        sampler = "saem"))))

  expect_s3_class(sampled, "ctJuliaFit")
  # And as a fit, not as a model spec: "ctFitModel" marks an unfitted handle,
  # and the optimised path gives its fits "ctFit".
  expect_s3_class(sampled, "ctFit")
  # The draws are where an optimised fit's normal-approximation draws live, so
  # everything downstream reads them without knowing which produced them.
  # Every draw of every chain, however many the effective-size target took:
  # 80 asked for is the size aimed at, not a fixed count (see the budget test
  # below), so the shape is checked against what was drawn.
  expect_equal(dim(sampled$estimate$rawposterior),
    c(2L * sampled$sample$draws, npar))
  expect_gte(sampled$sample$draws, 1L)
  expect_lte(sampled$sample$draws, 4L * 80L)
  expect_true(all(is.finite(sampled$estimate$rawposterior)))
  expect_equal(colnames(sampled$estimate$rawposterior),
    ctsem:::.ctBackendRawParameterNames(fit, npar))

  # The point estimate becomes the posterior mean, and the point the chains
  # were placed from is kept rather than overwritten: on the joint target
  # SAEM's estimate (it ran from the fit's), with what SAEM did.
  expect_equal(sampled$estimate$raw,
    as.numeric(colMeans(sampled$estimate$rawposterior)))
  expect_identical(sampled$sample$placement$method, "saem")
  expect_gt(sampled$sample$placement$saem_iterations, 0L)
  expect_length(sampled$estimate$placed_raw, npar)
  expect_true(all(is.finite(sampled$estimate$placed_raw)))
  expect_equal(sampled$sample$start, sampled$estimate$placed_raw)
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
  # intoverpop='none': this test asserts on the per-subject effect summaries
  # below, which the marginal target -- the default -- does not carry at all.
  # SAEM's kernel by name: its diagnostics (the scale moves' acceptance) are
  # among what is asserted, and on this fixture the joint default is the map.
  sampled <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 2, warmup = 80, draws = 80, intoverpop = "none",
        sampler = "saem"))))
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
  expect_length(diagnostics$accept, 2L * diagnostics$draws)
  # The point the chains started from.
  expect_equal(diagnostics$start, sampled$estimate$placed_raw)
  # SAEM's kernel, as asked.
  expect_identical(diagnostics$sampler, "saem")
  expect_length(diagnostics$scale_accept, 2L)
  expect_output(print(diagnostics), "SAEM kernel")

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

test_that("sampled effects are each parameter's natural deviation, at reduced rank and nested", {
  skip_without_julia()
  # What a sampled fit reports as its effects is each parameter's deviation
  # from its population value, per subject or group -- not the sampler's
  # standardised coordinates, which depend on a choice of square root of the
  # population covariance and, at reduced rank, are fewer than the parameters.
  data <- laplace_fixture_data(12L, 6L, 31L)
  data$study <- ceiling(data$id / 3)
  model <- laplace_fixture_model()
  model$pars$indvarying[model$pars$matrix == "T0MEANS" & is.na(model$pars$value)] <- TRUE
  warned <- character()
  # The joint route by name: 'auto' samples this model's exact marginal (its
  # effects shift means, with a Gaussian indicator), which has no effects to
  # save.
  run <- function(m, ...) withCallingHandlers(suppressMessages(ctFit(data, m,
    backend = "julia", cores = 1, priors = TRUE, optimize = FALSE,
    intoverpop = "none", ...,
    sampleControl = list(chains = 1L, warmup = 20L, draws = 20L, minESS = 0,
      saveEffects = TRUE, processes = FALSE))),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  # Rank one, a loading matrix over both effects -- on the joint route
  # ('none'), which builds it as laplace does: one coordinate per subject, two
  # parameters, so two deviations each.
  reduced <- run(model, poprank = 1)
  expect_length(reduced$sample$effect_mean, 12L * 2L)
  expect_false(is.null(names(reduced$sample$effect_mean)))
  expect_equal(dim(reduced$sample$effects), c(reduced$sample$draws, 24L))
  expect_true(all(is.finite(reduced$sample$effects)))
  # SAEM places these chains, and the optimisation before it runs no
  # iterations, so its verdict ("no progress ... Treat this fit as failed") is
  # kept on the fit rather than raised, and neither it nor that optimisation's
  # log likelihood is reported as the sampled fit's.
  expect_false(any(grepl("Treat this fit as failed", warned)))
  expect_true(any(grepl("no progress", reduced$optim$placement_warnings)))
  printed <- utils::capture.output(print(reduced))
  expect_false(any(grepl("^  converged:|log likelihood", printed)))
  s <- summary(reduced)
  expect_null(s$loglik)
  expect_null(s$aic)

  # Nested: a study-level effect on the manifest mean too -- 12 subjects x 2
  # plus 4 studies x 1, labelled from the engine's deviation layout.
  nested <- model
  nested$groupIDnames <- "study"
  nested$pars$indvarying <- ifelse(nested$pars$indvarying %in% TRUE,
    ifelse(nested$pars$matrix == "MANIFESTMEANS", TRUE, nested$subjectIDname),
    FALSE)
  twolevel <- run(nested)
  expect_length(twolevel$sample$effect_mean, 12L * 2L + 4L)
  expect_false(is.null(twolevel$sample$effectIndex))
  expect_equal(nrow(twolevel$sample$effectIndex), 28L)
  expect_identical(colnames(twolevel$sample$effects), twolevel$sample$effectIndex$label)
})

test_that("the effects come back summarised, or in full when asked for", {
  skip_without_julia()
  fit <- .sample_fixture()
  neffects <- length(fit$model_spec$subject_starts) *
    length(fit$model_spec$laplace$re_index)

  # intoverpop='none': there is no random effect in the sampled vector at all
  # under the marginal default, so `saveEffects` has nothing to save -- this
  # test is specifically about that vector.
  #
  # By default only the summary crosses the bridge: the draws themselves are
  # nsubjects x neffects x chains x draws numbers, and the bridge moves about
  # 1 MB/s.
  lean <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 1, warmup = 60, draws = 60, intoverpop = "none"))))
  expect_null(lean$sample$effects)
  expect_length(lean$sample$effect_mean, neffects)
  expect_length(lean$sample$effect_sd, neffects)
  expect_true(all(lean$sample$effect_sd > 0))

  full <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 1, warmup = 60, draws = 60, saveEffects = TRUE,
        intoverpop = "none"))))
  expect_equal(dim(full$sample$effects), c(full$sample$draws, neffects))
  expect_true(all(is.finite(full$sample$effects)))
})

test_that("NUTS samples the joint posterior through the same runner", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  neffects <- length(fit$model_spec$subject_starts) *
    length(fit$model_spec$laplace$re_index)
  # control$sampler = 'nuts': NUTS on the parameters and every effect at once,
  # the alternative to the default SAEM kernel on the joint posterior
  # (`ctsem_sample`, checked against SAEM's kernel in the engine suite),
  # placed and stopped as the default is. What is left for here is that it
  # reaches the same assembly with the same shape.
  nuts <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 2, warmup = 80, draws = 80, intoverpop = "none",
        sampler = "nuts", processes = FALSE))))
  d <- nuts$sample
  expect_identical(d$sampler, "nuts")
  expect_identical(d$target, "joint")
  expect_equal(dim(nuts$estimate$rawposterior), c(2L * d$draws, npar))
  expect_true(all(is.finite(nuts$estimate$rawposterior)))
  expect_length(d$rhat, npar)
  expect_length(d$stepsize, 2L)
  expect_true(all(is.finite(d$ebfmi)))
  expect_length(d$accept, 2L * d$draws)
  expect_length(d$effect_mean, neffects)
  expect_true(all(d$effect_sd > 0))
  expect_output(print(d), "ctsem Hamiltonian sample")
  # Placed by SAEM's state by default on the joint target, whichever sampler.
  expect_identical(d$placement$method, "saem")
  expect_equal(d$start, nuts$estimate$placed_raw)
})

# The map sampler on one fit, and what it is checked against. Its change of
# variables first, deterministically (`ctsem_sample_map_check`: the latent
# density's gradient and the map's log Jacobian against differences), then its
# draws: the parameters' means against `against`'s, as z on both runs'
# effective sizes, and the reported effects' means against SAEM's kernel,
# which samples every effect. Widths are not compared: at these sizes the
# other samplers mix slowly along this fixture's ridge between the initial and
# manifest means (effective sizes 16-73 there, against the map's 215-634, and
# SAEM's sds 0.065 against the exact marginal's 0.138), so a width ratio reads
# the reference's error, not the map's.
# The laplace fixture's initial and manifest means trade off along a ridge that
# only the trajectories' decay separates: their sum is identified, their split
# barely. At four chains of 1,500 draws SAEM and the map each put the split's
# sd anywhere from 0.065 to 0.23 from run to run, with effective sizes of
# hundreds that do not show it (`scripts/marginal-default-2026-10/
# map_start_diag.R`), and the sum the same in every run. So the means are
# compared with the pair's sum in place of the pair.
.map_ridge <- function(X) {
  pair <- c("T0m_eta1", "mm_Y1")
  if (!all(pair %in% colnames(X))) return(X)
  cbind(X[, setdiff(colnames(X), pair), drop = FALSE],
    ridge_sum = X[, pair[1]] + X[, pair[2]])
}

.map_checked <- function(fit, against, draws = 300L) {
  ask <- function(...) suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1, control = list(chains = 2, draws = draws,
      processes = FALSE, ...))))
  map <- ask(intoverpop = "none", sampler = "map", saveEffects = TRUE)
  ref <- do.call(ask, against)
  saem <- ask(intoverpop = "none")
  a <- .map_ridge(map$estimate$rawposterior)
  b <- .map_ridge(ref$estimate$rawposterior)
  sb <- apply(b, 2, stats::sd)
  # The pair's own effective sizes, the smaller, for their sum.
  ess <- function(r) {
    e <- pmax(as.numeric(r$sample$ess), 1)
    names(e) <- colnames(r$estimate$rawposterior)
    keep <- e[setdiff(names(e), c("T0m_eta1", "mm_Y1"))]
    if (length(keep) == length(e)) e else
      c(keep, ridge_sum = min(e[c("T0m_eta1", "mm_Y1")]))
  }
  check <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaModule(fit$model_spec$project)$
    ctsem_sample_map_check(ctsem:::.ctJuliaObjective(fit, hybrid = TRUE),
      as.numeric(fit$estimate$raw),
      hessian = ctsem:::.ctJuliaPut(as.matrix(fit$uncertainty$hessian))))
  list(map = map, saem = saem, check = check,
    z = max(abs(colMeans(a) - colMeans(b)) / (sb * sqrt(1 / ess(map) + 1 / ess(ref)))),
    effects = max(abs(map$sample$effect_mean - saem$sample$effect_mean) /
      saem$sample$effect_sd))
}

test_that("the map sampler draws the joint posterior with the filter's effects integrated", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  neffects <- length(fit$model_spec$subject_starts) *
    length(fit$model_spec$laplace$re_index)
  # control$sampler = 'map' (`ctsem_sample_map`): NUTS through a learnt map on
  # the hybrid objective. The fixture's one effect is on a manifest mean, which
  # the filter integrates exactly, so no unit keeps a coordinate: the chains
  # move the parameters alone, and every effect is reported from its exact
  # conditional given each draw. Every unit being quadratic, the Laplace
  # marginal is exact, and NUTS on it is the reference for the parameters.
  r <- .map_checked(fit, list(intoverpop = "laplace", sampler = "nuts"))
  expect_equal(r$check$coordinates, 0L)
  expect_lt(r$check$gradient, 1e-5)
  expect_lt(r$check$jacobian, 1e-6)
  d <- r$map$sample
  expect_identical(d$sampler, "map")
  expect_identical(d$target, "joint")
  expect_identical(d$placement$method, "fit")
  expect_equal(dim(r$map$estimate$rawposterior), c(2L * d$draws, npar))
  expect_true(all(is.finite(r$map$estimate$rawposterior)))
  expect_length(d$effect_mean, neffects)
  expect_identical(names(d$effect_mean), names(r$saem$sample$effect_mean))
  expect_equal(ncol(d$effects), neffects)
  expect_true(all(is.finite(d$effects)))
  expect_true(all(d$effect_sd > 0))
  # A wrong Jacobian or conditional moves a mean by many Monte Carlo errors
  # and an effect's mean by a good part of its sd; these see that.
  expect_lt(r$z, 4)
  expect_lt(r$effects, 0.5)
  expect_match(summary(r$map)$uncertaintyNote, "NUTS through a learnt map", fixed = TRUE)

  # Refused where it cannot apply, by name.
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "map", intoverpop = "laplace")),
    "control$sampler = 'map'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "map", intoverpop = "none", placement = "saem")),
    "does not apply to control$sampler = 'map'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "map", intoverpop = "none", adapt_effects = TRUE)),
    "control$adapt_effects", fixed = TRUE)
})

test_that("the map sampler maps the effects the filter cannot integrate", {
  skip_without_julia()
  # A random drift beside the manifest mean: the drift is not affine in the
  # filter, so it stays a coordinate, partially non-centred through the
  # subjects' summaries -- and the manifest mean with it, one qualifying
  # effect beside another costing more in the filter than as a coordinate
  # (the hybrid's cost rule). No exact marginal here, so the parameters'
  # means are read against SAEM's.
  fit <- fit_cached("map fixture: random drift", {
    model <- laplace_fixture_model()
    model$pars$indvarying[model$pars$matrix == "DRIFT"] <- TRUE
    suppressWarnings(suppressMessages(ctsem::ctFit(laplace_fixture_data(), model,
      backend = "julia", cores = 1, intoverpop = "laplace", priors = TRUE,
      optimcontrol = list(finishsamples = 20))))
  })
  r <- .map_checked(fit, list(intoverpop = "none"))
  expect_gte(r$check$coordinates, length(fit$model_spec$subject_starts))
  expect_lt(r$check$gradient, 1e-5)
  expect_lt(r$check$jacobian, 1e-6)
  expect_identical(r$map$sample$sampler, "map")
  expect_length(r$map$sample$effect_mean, length(r$saem$sample$effect_mean))
  expect_lt(r$z, 4)
  expect_lt(r$effects, 0.5)
})

test_that("placement = 'fit' starts the chains around the fit's own estimate", {
  skip_without_julia()
  fit <- .sample_fixture()
  # The old placement, kept as a choice: no SAEM run, the chains drawn around
  # the fit's estimate.
  sampled <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 1, warmup = 30, draws = 30, placement = "fit"))))
  expect_identical(sampled$sample$placement$method, "fit")
  expect_equal(sampled$estimate$placed_raw, as.numeric(fit$estimate$raw))
  # A marginal target has no effects for SAEM to start, so the default there
  # is 'fit', and asking for 'saem' is refused by name.
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(intoverpop = "laplace", placement = "saem")),
    "control$placement = 'saem'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(placement = "laplace")), "'saem' or 'fit'", fixed = TRUE)
})

test_that("the SAEM kernel is refused by name where it cannot apply", {
  skip_without_julia()
  fit <- .sample_fixture()
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "saem", intoverpop = "laplace")),
    "control$sampler = 'saem'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "saem", intoverpop = "none", adapt_metric = TRUE)),
    "control$adapt_metric", fixed = TRUE)
  # Refused under the default too, since that is SAEM's kernel on the joint
  # posterior, and the message says how to get NUTS.
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(intoverpop = "none", settleTol = 0.1)),
    "sampler = 'nuts'", fixed = TRUE)
  # The default follows the target: the marginal one has no effects for
  # SAEM's sweeps.
  expect_identical(.ctBackendSamplerName(NULL, joint = TRUE), "saem")
  # The joint default is the map where the filter integrates some of the
  # effects, as it does this fixture's one.
  expect_identical(.ctBackendSamplerName(NULL, joint = TRUE, integrable = TRUE), "map")
  expect_true(.ctBackendMapDefault(.ctBackendSpec(fit)))
  expect_false(.ctBackendMapDefault(.ctBackendSpec(fit), tisampled = TRUE))
  expect_identical(.ctBackendSamplerName(NULL, joint = FALSE), "nuts")
  expect_identical(.ctBackendSamplerName("nuts", joint = TRUE), "nuts")
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "hmc")), "'nuts', 'saem', 'map' or 'smc'", fixed = TRUE)
})

test_that("a marginal target's default sampler is chosen by the fit's check", {
  # Tempered SMC where the check named a direction in which the posterior
  # departs from the Hessian's normal, NUTS where it named none or where the
  # call set a control only NUTS reads; the joint target keeps SAEM's kernel.
  # The tests below of NUTS's own mechanics on a marginal target -- its
  # worker processes, step size, draw budget, R-hat -- name it for that
  # reason: their fixtures' checks name directions.
  name <- ctsem:::.ctBackendSamplerName
  expect_identical(name(NULL), "saem")
  expect_identical(name(NULL, joint = TRUE, named = 3L), "saem")
  expect_identical(name(NULL, joint = FALSE), "nuts")
  expect_identical(name(NULL, joint = FALSE, named = 3L), "smc")
  expect_identical(name(NULL, joint = FALSE, named = 3L, nutsonly = TRUE), "nuts")
  expect_identical(name("nuts", joint = FALSE, named = 3L), "nuts")
  expect_identical(name("smc", joint = FALSE), "smc")
  # What counts as asking for NUTS: its chain count, warmup or starting step,
  # or a setting of its warmup metric; not the draws, which SMC reads too.
  nutsctl <- ctsem:::.ctBackendNutsControls
  asked <- ctsem:::.ctBackendNutsAsked
  expect_identical(nutsctl(list()), character(0))
  expect_identical(nutsctl(list(draws = 10, chains = 2, seed = 1)), "chains")
  expect_identical(nutsctl(list(stepsize = 0.1, warmup = 50, chains = NULL)),
    c("warmup", "stepsize"))
  expect_false(asked(list(draws = 400, minESS = 100, particles = 40)))
  expect_true(asked(list(chains = 2)))
  expect_true(asked(list(adapt_metric = TRUE)))
  # ctFit(optimize = FALSE) folds its own resolved chains into every control,
  # and says which the caller set instead.
  expect_false(asked(list(chains = 4, warmup = 200), nutsgiven = character(0)))
})

test_that("tempered SMC samples the marginal posterior through the same assembly", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  # control$sampler = 'smc' (`ctsem_sample_smc`, checked on closed forms in
  # the engine suite): particles carried from the check's normal to the
  # posterior, then as many of them chains as the effective size needs. What
  # is left for here is the assembly, the
  # record, the refusals, and agreement with NUTS on the same marginal.
  smc <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(intoverpop = "laplace", sampler = "smc", particles = 40,
        draws = 400, processes = FALSE))))
  d <- smc$sample
  expect_identical(d$sampler, "smc")
  expect_identical(d$target, "marginal")
  expect_identical(d$nparticles, 40L)
  expect_true(d$chains >= 16L && d$chains <= 40L)
  expect_gte(d$final_tau, 1)
  expect_identical(d$warmup, 0L)
  expect_gte(d$draws, 20L)
  expect_equal(dim(smc$estimate$rawposterior), c(d$chains * d$draws, npar))
  expect_true(all(is.finite(smc$estimate$rawposterior)))
  expect_gte(d$nstages, 1L)
  expect_equal(utils::tail(d$stage_beta, 1), 1)
  expect_length(d$stepsize, d$chains)
  expect_true(all(d$stepsize == d$stepsize[1]))
  expect_true(is.infinite(d$rhat_target))
  # No R-hat: its chains start from particles already drawn from the
  # posterior, so their agreement says nothing about mixing.
  expect_length(d$rhat, npar)
  expect_true(all(is.na(d$rhat)))
  expect_true(all(is.na(summary(smc)$popmeans$Rhat)))
  expect_output(print(d), "tempered SMC")
  # The other sampler on the same target, for the means.
  nuts <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(intoverpop = "laplace", chains = 2, warmup = 100, draws = 400,
        processes = FALSE))))
  expect_lt(max(abs(colMeans(smc$estimate$rawposterior) -
    colMeans(nuts$estimate$rawposterior)) / nuts$estimate$se), 0.5)
  # Refused by name on the joint target, which has no version of it.
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "smc", intoverpop = "none")),
    "control$sampler = 'smc'", fixed = TRUE)
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "smc", intoverpop = "laplace", settleTol = 0.1)),
    "sampler = 'nuts'", fixed = TRUE)
  # As are NUTS's own counts beside it, which SMC would otherwise ignore.
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(sampler = "smc", intoverpop = "laplace", chains = 2)),
    "control$chains is a setting of NUTS", fixed = TRUE)
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
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 1, warmup = 40, draws = 40))),
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

test_that("the marginal sampler starts from the normal check, made where the fit has none", {
  skip_without_julia()
  fit <- .sample_fixture()
  expect_false(is.null(fit$uncertainty$details$normalcheck))
  # A fit that has its check keeps it, gaining only DAIS's normal where the
  # check named directions and the fit has no posterior covariance of its own;
  # one without gets 'checked''s; one whose climb took a jump in keeps its
  # Hessian.
  kept <- ctsem:::.ctBackendSampleCheck(fit)
  expect_identical(kept$uncertainty$details$normalcheck,
    fit$uncertainty$details$normalcheck)
  added <- setdiff(names(kept$uncertainty$details), names(fit$uncertainty$details))
  expect_true(all(added %in% "dais"))
  # A default fit carries the curvature's normal and the check beside it, and
  # where the check names a direction DAIS is made for the metric, as
  # 'checked' would have made it.
  named <- fit
  named$uncertainty$details[c("importance_sampling", "dais")] <- NULL
  named$uncertainty$details$normalcheck$directions <-
    diag(length(fit$estimate$raw))[, 1L, drop = FALSE]
  named$uncertainty$details$normalcheck$widths <- 1.5
  set.seed(2)
  daised <- ctsem:::.ctBackendSampleCheck(named)
  expect_identical(ctsem:::.ctBackendMetricSource(daised), "dais")
  bare <- fit
  bare$uncertainty$details$normalcheck <- NULL
  set.seed(1)
  made <- ctsem:::.ctBackendSampleCheck(bare)$uncertainty$details$normalcheck
  expect_true(all(c("directions", "widths", "centre") %in% names(made)))
  expect_gt(made$draws, 0)
  jumped <- bare
  jumped$laplace$jump <- list(units = 1L)
  expect_identical(ctsem:::.ctBackendSampleCheck(jumped), jumped)

  # The sampler asks for it on the marginal target, not on the joint one,
  # whose engine builds its own population block.
  asked <- list()
  testthat::local_mocked_bindings(
    .ctBackendSampleCheck = function(fit, cores = 1L, verbose = 0) {
      asked[[length(asked) + 1L]] <<- cores
      fit
    },
    .ctBackendSampleRun = function(fit, target, ...) stop("no sampling here"))
  for (route in c("laplace", "none")) expect_error(ctFitUncertainty(bare, "sample",
    cores = 1, control = list(intoverpop = route)), "no sampling here")
  expect_length(asked, 1L)
})

test_that("NUTS on a marginal target records the coreset its chains moved on", {
  skip_without_julia()
  fit <- .sample_fixture()
  # Twelve units: each chain fits a coreset from its warm-up where one passes,
  # and its draws' ends are tested against the full density
  # (sample_surrogate.jl, checked against closed forms in the engine suite).
  # Here, that the record says per chain what happened, and that the oracle
  # records nothing.
  run <- function(coreset) suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(sampler = "nuts", chains = 2, warmup = 60, draws = 60,
        processes = FALSE, minESS = 0, coreset = coreset))))
  core <- run(TRUE)
  s <- core$sample
  expect_length(s$chain_coreset_size, 2L)
  expect_identical(s$coreset, any(s$chain_coreset_size > 0))
  expect_true(all(s$chain_coreset_accepted <= s$chain_coreset_proposed))
  expect_true(all(s$chain_coreset_size == 0 | s$chain_coreset_size <= s$nunits))
  expect_output(print(s), "ctsem Hamiltonian sample")
  expect_true(all(is.finite(core$estimate$rawposterior)))
  full <- run(FALSE)
  expect_false(full$sample$coreset)
  expect_true(all(full$sample$chain_coreset_size == 0))
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
  #
  # `placement = 'fit'`: under the default SAEM placement each worker runs its
  # own SAEM for its block, so the two routes start from different points
  # and cannot reproduce each other. What this checks -- the pooling and the
  # streams -- is the same under either placement.
  #
  # Four chains on two workers, two to a block, so a block's second chain has
  # to find its stream from the block's seed. `minESS = 0`: every chain draws
  # exactly the count asked for, so the two routes draw the same chains whole.
  run <- function(processes, minESS) suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 2,
      control = list(chains = 4, warmup = 0, draws = 3, seed = 777,
        processes = processes, placement = "fit", minESS = minESS))))
  inprocess <- run(FALSE, 0)
  viaprocess <- run(TRUE, 0)

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

  expect_identical(viaprocess$sample$draws, inprocess$sample$draws)
  expect_equal(dim(viaprocess$estimate$rawposterior),
    dim(inprocess$estimate$rawposterior))
  expect_equal(viaprocess$estimate$rawposterior,
    inprocess$estimate$rawposterior, tolerance = 1e-6)

  # With the default effective-size target, which three draws cannot meet, the
  # run is extended to its budget: in one session every chain by the same
  # amount, across the workers by amounts in proportion to each worker's
  # speed, so the chains may end at different lengths. Every draw is kept --
  # the pool once cut every chain to the shortest -- and each chain's first
  # batch is still the in-session chain's.
  targeted <- run(TRUE, 200)
  inprocessTargeted <- run(FALSE, 200)
  expect_gt(inprocessTargeted$sample$draws, 3L)
  skip_if_not(isTRUE(targeted$sample$processes),
    "processes = TRUE fell back to in-process sampling in this session")
  lengths <- ctsem:::.ctSampleChainDraws(targeted$sample)
  expect_length(lengths, 4L)
  expect_identical(sum(lengths), nrow(targeted$estimate$rawposterior))
  expect_true(all(lengths > 3L))
  first <- function(f) {
    L <- ctsem:::.ctSampleChainDraws(f$sample)
    starts <- cumsum(c(0L, L[-length(L)]))
    f$estimate$rawposterior[as.vector(outer(1:3, starts, "+")), , drop = FALSE]
  }
  expect_equal(first(targeted), first(inprocessTargeted), tolerance = 1e-6)
  expect_length(targeted$sample$rhat, ncol(targeted$estimate$rawposterior))
  expect_true(all(is.finite(targeted$sample$ess)))
})

test_that("intoverpop alone says what is sampled", {
  route <- ctsem:::.ctBackendSampleRoute
  # A route the caller named is sampled as named; a Laplace route 'auto'
  # chose to optimise is sampled as 'none', the joint posterior, as
  # ctFit(optimize = FALSE, intoverpop = 'auto') samples it. Before
  # 2026-10-07 a separate `target` said the same thing a second way.
  expect_identical(route("laplace", requested = "laplace"), "laplace")
  expect_identical(route("laplace", requested = TRUE), "laplace")
  expect_identical(route("laplace", requested = "auto"), "none")
  expect_identical(route("laplace", requested = NULL), "none")
  # Unless the Laplace term is exact (every unit's integrand Gaussian in its
  # effects): then its marginal is the exact posterior's, and that is sampled.
  expect_identical(route("laplace", requested = "auto", exact = TRUE), "laplace")
  expect_identical(route("laplace", requested = NULL, exact = TRUE), "laplace")
  expect_identical(route("none", requested = "auto", exact = TRUE), "none")
  expect_identical(route("laplace", requested = "auto", override = "none",
    exact = TRUE), "none")
  expect_identical(route("augmented", requested = "auto"), "augmented")
  expect_identical(route("none", requested = FALSE), "none")
  # A backend caller passing only the route has named it.
  expect_identical(route("laplace"), "laplace")
  # ctFitUncertainty's control$intoverpop overrides, in the same words, within
  # what the fit's specification can carry.
  expect_identical(route("laplace", requested = "laplace", override = "none"), "none")
  expect_identical(route("laplace", requested = "auto", override = "laplace"), "laplace")
  expect_identical(route("laplace", requested = "laplace", override = "auto"), "laplace")
  expect_error(route("augmented", override = "none"), "Refit with intoverpop = 'none'",
    fixed = TRUE)
  expect_error(route("laplace", override = "augmented"), "needs a fit made with",
    fixed = TRUE)
  expect_error(route("laplace", override = "joint"), "'auto', 'none', 'laplace'",
    fixed = TRUE)
  # The retired name, and ctFit's own argument given a second time, are refused
  # by name rather than accepted beside intoverpop.
  expect_error(ctsem:::.ctBackendSampleCheckControl(list(target = "joint")),
    "intoverpop says what is sampled", fixed = TRUE)
  expect_error(ctsem:::.ctBackendSampleCheckControl(list(intoverpop = "none"),
    listname = "sampleControl"), "ctFit()'s own intoverpop", fixed = TRUE)
  expect_null(ctsem:::.ctBackendSampleCheckControl(list(intoverpop = "none")))
})

test_that("ctFit(optimize = FALSE) samples the exact marginal where the Laplace term is exact", {
  skip_without_julia()
  data <- laplace_fixture_data()
  prep <- function(model, ...) suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, optimize = FALSE, fit = FALSE, ...)))
  route <- function(x) ctsem:::.ctBackendIntOverPop(ctsem:::.ctBackendSpec(x))
  # A random manifest mean with a Gaussian indicator: the Laplace term is the
  # exact integral, so 'auto' samples that marginal, and the fit says why.
  exact <- prep(laplace_fixture_model())
  expect_identical(route(exact), "laplace")
  expect_identical(exact$args$resolved$intoverpop, "laplace")
  expect_match(exact$args$resolved$intoverpopreason, "exact marginal", fixed = TRUE)
  # Named, the joint is the joint.
  expect_identical(route(prep(laplace_fixture_model(), intoverpop = "none")), "none")
  # A random drift is integrated by an approximation, so 'auto' samples the
  # joint posterior.
  drift <- laplace_fixture_model()
  drift$pars$indvarying <- drift$pars$matrix == "DRIFT"
  approx <- prep(drift)
  expect_identical(route(approx), "none")
  expect_identical(approx$args$resolved$intoverpop, "none")
})

test_that("ctFit(optimize = FALSE, intoverpop = 'laplace') samples the Laplace marginal in worker processes", {
  skip_without_julia()
  # Each worker builds the Laplace objective for its own chain; that the
  # pooled draws reproduce the in-session ones on this target is the test
  # above, whose fixture names 'laplace'. This is the ctFit() entry, which
  # sampled the joint posterior until 2026-10-07.
  fit <- suppressWarnings(suppressMessages(ctFit(laplace_fixture_data(),
    laplace_fixture_model(), backend = "julia", cores = 2,
    intoverpop = "laplace", priors = TRUE, optimize = FALSE,
    sampleControl = list(sampler = "nuts", chains = 2, warmup = 30, draws = 30, minESS = 0,
      processes = TRUE))))
  expect_identical(fit$sample$target, "marginal")
  expect_identical(fit$sample$sampler, "nuts")
  expect_identical(fit$args$resolved$intoverpop, "laplace")
  expect_identical(ncol(fit$estimate$rawposterior), length(fit$estimate$raw))
  expect_true(all(is.finite(fit$estimate$rawposterior)))
  skip_if_not(isTRUE(fit$sample$processes),
    "processes = TRUE fell back to in-process sampling in this session")
  # A chain from each worker, each to the count asked for.
  expect_identical(nrow(fit$estimate$rawposterior), 60L)
})

test_that("sampling workers start Julia at their share of the cores", {
  skip_without_julia()
  skip_if_not(ctsem:::.ctBackendCanWarm(), "future is not installed")
  fit <- .sample_fixture()
  # A Julia session's thread count is fixed when it starts. Workers were warmed
  # at one thread whatever a chain later asked for, so processes mode at 16
  # cores and four chains ran on four (bigre: 352 s against 106 s for the same
  # chains as threads in one session, dev1), and nothing noticed.
  nthreads <- function(handles) {
    skip_if(ctsem:::.ctBackendWarmWait(handles) < 1L,
      "no worker could be warmed in this session")
    future::value(future::future(as.integer(utils::getFromNamespace(
      ".ctJuliaEval", "ctsem")("Threads.nthreads()")), seed = TRUE))
  }
  on.exit(ctsem:::.ctBackendWarmStop(NULL), add = TRUE)
  expect_equal(nthreads(ctsem:::.ctBackendWarmWorkers(fit, workers = 2L,
    threads = 2L)), 2L)
  # A pool warmed at another width is replaced, not reused.
  expect_equal(nthreads(ctsem:::.ctBackendWarmWorkers(fit, workers = 2L,
    threads = 1L)), 1L)
})

test_that("chains in worker processes stop together, by the run's rule", {
  skip_without_julia()
  fit <- .sample_fixture()
  npar <- length(fit$estimate$raw)
  # Workers that stopped on their own share of the target returned chains of
  # different lengths, and the pool cut them all to the shortest: 230 of 930,
  # 1365, 230 and 696 draws on one run (dev2, 2026-10-07), below the target.
  # Now a worker judges nothing and the parent decides every batch.
  wc <- ctsem:::.ctBackendWorkerControl(list(minESS = 400, meanESS = 800,
    maxDraws = 50L, chains = 4L))
  expect_identical(wc$minESS, 0)
  expect_null(wc$meanESS)
  expect_null(wc$maxDraws)
  expect_identical(wc$chains, 4L)
  # A run without a target has nothing to decide.
  expect_null(ctsem:::.ctBackendChainCoordinator(fit,
    ctsem:::.ctBackendSampleControl(list(minESS = 0)), list(first = 20L), 4L, 2L))

  settings <- ctsem:::.ctBackendSampleControl(list(minESS = 1000))
  budget <- ctsem:::.ctBackendSampleBudget(100L, 4L, settings)
  first <- as.integer(budget$first)
  co <- ctsem:::.ctBackendChainCoordinator(fit, settings, budget, chains = 4L,
    workers = 2L)
  on.exit(co$close(), add = TRUE)
  set.seed(3)
  # A worker's batch as `.ctBackendChainCoordinate()` writes it: two chains
  # each, the first worker's running twice as fast as the second's.
  batch <- function(rate = 1) list(draws = matrix(stats::rnorm(npar * 2L * first), npar),
    total = first, nchains = 2L, rate = rate)
  verdict <- function(r, w = 1L) file.path(co$dir, sprintf("verdict_w%d_r%d.rds", w, r))
  # One worker of two has written its batch: nothing is decided yet.
  # Each tick is called outside the expectation: it acts, and an expectation
  # may evaluate its argument more than once.
  saveRDS(batch(2), file.path(co$dir, "w1_r1.rds"))
  said <- co$tick(c(FALSE, FALSE))
  expect_null(said)
  expect_false(file.exists(verdict(1L)))
  # Both have: every chain is judged at once, and every chain is asked for the
  # same number more -- 400 independent draws are well short of ESS 1000.
  saveRDS(batch(), file.path(co$dir, "w2_r1.rds"))
  said <- co$tick(c(FALSE, FALSE))
  expect_match(said, paste(first, "draws per chain"), fixed = TRUE)
  wanted <- readRDS(verdict(1L))
  expect_gt(wanted, 0L)
  # Shared by speed: the faster worker's chains are asked for twice as many,
  # so both finish the round together.
  expect_equal(wanted / readRDS(verdict(1L, 2L)), 2, tolerance = 0.05)
  expect_identical(ctsem:::.ctBackendShareDraws(600L, c(2, 1), c(2L, 2L)),
    c(200L, 100L))
  expect_identical(ctsem:::.ctBackendShareDraws(600L, c(NA, 1), c(2L, 2L)),
    c(150L, 150L))
  expect_identical(ctsem:::.ctBackendWorkerThreads(6L, 4L), c(2L, 2L, 1L, 1L))
  expect_identical(ctsem:::.ctBackendWorkerThreads(16L, 4L), rep(4L, 4L))
  expect_identical(ctsem:::.ctBackendWorkerThreads(2L, 4L), rep(1L, 4L))
  # A worker's half writes its batch and returns the parent's answer.
  answered <- ctsem:::.ctBackendChainCoordinate(co$dir, 1L)(batch()$draws, first)
  expect_identical(answered, wanted)
  # A worker that ended without writing stops the others rather than leaving
  # them waiting on a batch that is never coming.
  saveRDS(batch(), file.path(co$dir, "w1_r2.rds"))
  said <- co$tick(c(FALSE, TRUE))
  expect_null(said)
  expect_identical(readRDS(verdict(2L)), 0L)
  expect_identical(readRDS(verdict(2L, 2L)), 0L)
  # And a worker still waiting when the parent finishes is released.
  co$close()
  answered <- ctsem:::.ctBackendChainCoordinate(co$dir, 2L)(batch()$draws, first)
  expect_identical(answered, 0L)
})

test_that("chains sampled in this session run one after another, whatever the pool holds", {
  skip_without_julia()
  fit <- .sample_fixture()
  # With one thread there is no second worker to put a chain on, and nothing
  # here can fail.
  skip_if(!isTRUE(ctsem:::.ctBackendSessionThreads() >= 2L),
    "the julia session has a single thread")

  # An evaluation outside any `cores` ceiling spreads over every thread and
  # leaves the worker pool at full width, its workers idle -- where a fit's
  # random-effect check once left it too. In-process chains share one
  # objective, so a chain runner that handed the second chain to an idle
  # worker had two chains writing the same workspaces at once, and the draws
  # came out order-1 different from the same seeds.
  invisible(ctJuliaEvaluate(fit))
  # placement = 'fit': SAEM's placement draws every chain's start from one
  # stream seeded by the run, so under it chain 2 of seed s is not the chain
  # of seed s + 1 and this identity cannot hold. Around the fit, each chain's
  # stream is its own.
  draws <- function(chains, seed) suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 2,
      control = list(chains = chains, warmup = 0, draws = 3, seed = seed,
        processes = FALSE, placement = "fit"))))$estimate$rawposterior
  both <- draws(2L, 777L)
  # Chain `c` of a run seeded `s` draws the stream `s + c`, so the second
  # chain of the pair is the only chain of a run seeded 778.
  expect_equal(both, rbind(draws(1L, 777L), draws(1L, 778L)),
    tolerance = 1e-6)
})

test_that("ctFitUncertainty(fit, 'sample') refuses what it cannot sample", {
  skip_without_julia()
  expect_error(ctFitUncertainty(list(), uncertainty = "sample"),
    "backend='julia'")

  # A stan fit is refused by name, pointing at its own sampling entry point
  # rather than at a julia-only error about model structure.
  expect_error(ctFitUncertainty(ctstantestfit, uncertainty = "sample"),
    "backend='stan'.*optimize.*FALSE")

  # The augmented route has no separate posterior over the random effects --
  # they are carried in the state, integrated by the filter -- so the joint
  # target is refused by name. The marginal is not: it is the same augmented
  # marginal `ctFit(optimize=FALSE, intoverpop='augmented')` already samples, and
  # `uncertainty = 'sample'` reaches it by the same route.
  set.seed(4)
  data <- do.call(rbind, lapply(1:8, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * 0.5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  augmented <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, optimcontrol = list(estonly = TRUE),
    intoverpop = "augmented")))
  expect_error(ctFitUncertainty(augmented, uncertainty = "sample",
    control = list(intoverpop = "none")), "made with intoverpop = 'augmented'",
    fixed = TRUE)

  marginal <- suppressWarnings(suppressMessages(ctFitUncertainty(augmented,
    uncertainty = "sample", cores = 1,
    control = list(chains = 1, warmup = 20, draws = 20))))
  expect_s3_class(marginal, "ctJuliaFit")
  expect_false(is.null(marginal$sample))
  expect_null(marginal$sample$effects)
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
  #
  # Called by its 3.11.1 name deliberately, since this test is also what
  # confirms that name still means exactly what it did once it became a thin
  # wrapper over ctFitUncertainty().
  fit <- .sample_fixture()
  sampled <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 2, warmup = 40, draws = 40))))
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

  # The Laplace fit this was sampled from is a genuinely optimized fit
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
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(sampler = "nuts", chains = 2, warmup = 120, draws = 120))))
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
  expect_match(summarised$sampleNote,
    paste0("2 chains x ", sampled$sample$draws, " draws"), fixed = TRUE)
  expect_identical(names(summarised)[1L], "sampleNote")
  # A sampled fit did not run an uncertainty pass, and used to say it had.
  # Naming the kernel that drew them: NUTS on this fixture's marginal target.
  expect_match(summarised$uncertaintyNote, paste0("posterior draws (",
    if (identical(sampled$sample$sampler, "saem")) "SAEM kernel" else "NUTS",
    ")"), fixed = TRUE)
  expect_false(grepl("ctFitUncertainty", summarised$uncertaintyNote, fixed = TRUE))

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
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(sampler = "nuts", chains = 2, warmup = 12, draws = 12))))
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

test_that("the process-path poll stops waiting on a worker that has gone", {
  skip_if_not_installed("future")
  # `resolved()` erroring is one way `future` reports a worker process that has
  # died. The poll has to count that chain as over and let the `value()` pass
  # after it report the failure, rather than stop the whole sample -- and
  # rather than read it as "not yet", which waits forever.
  testthat::local_mocked_bindings(
    resolved = function(x, ...) stop("simulated: the worker process has gone"),
    .package = "future")
  gone <- structure(list(), class = "Future")
  expect_null(ctsem:::.ctBackendReportProcesses(list(gone, gone),
    tempfile(c("chain1_", "chain2_")), chains = 2L, interval = 0.01,
    overwrite = FALSE))
})

test_that("an effective-size target turns the draw count into a budget", {
  skip_without_julia()
  fit <- .sample_fixture()

  # `minESS` used to do nothing at all without `maxDraws`: the engine extends
  # towards a budget that defaulted to exactly the draws asked for, so there
  # was nothing to extend into and the target could only be reported after the
  # fact. Now a target met early stops the run, and one not met extends it to
  # four times the count asked for (`.ctBackendSampleBudget()`).
  #
  # `rhatTarget` is raised out of the way because this fixture is twelve
  # subjects and will not reach 1.01; what is under test is the stopping rule,
  # not whether a small fixture mixes.
  met <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1,
    control = list(sampler = "nuts", chains = 2, warmup = 100, draws = 200, processes = FALSE,
      minESS = 0.5, rhatTarget = 100))))
  expect_lt(met$sample$draws, 200L)
  expect_gt(met$sample$draws, 0L)

  # `minESS = 0` is the off switch, and then it takes exactly what it was
  # asked for. Needed explicitly now that 200 is the default target -- a
  # default fit stops early, which is the point of it.
  full <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1,
    control = list(sampler = "nuts", chains = 2, warmup = 100, draws = 200, processes = FALSE,
      minESS = 0))))
  expect_equal(full$sample$draws, 200L)

  # And the default is a target with a budget of four times the count asked
  # for: the twelve-subject fixture does not reach min ESS 200 in 200 draws,
  # so the run goes on past them, stops at the budget at most, and says it
  # fell short.
  # Warned about as well; the recorded diagnosis is what is asserted.
  default <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1,
    control = list(sampler = "nuts", chains = 2, warmup = 100, draws = 200, processes = FALSE))))
  expect_gt(default$sample$draws, 200L)
  expect_lte(default$sample$draws, 800L)
  expect_identical(default$sample$ess_target, 200)
  if (min(default$sample$ess) < 200) {
    expect_match(paste(default$sample$diagnosis, collapse = "; "),
      "(target 200)", fixed = TRUE)
  }

  # The old spelling is refused with the new one named, rather than dropped by
  # `$` and silently ignored.
  expect_error(ctFitUncertainty(fit, uncertainty = "sample",
    control = list(minEss = 100)), "did you mean minESS")
})

test_that("a fixed stepsize is what every chain starts from", {
  skip_without_julia()
  fit <- .sample_fixture()

  # With no warmup there is nothing to move the step size from where it began,
  # so this is the one setting where the initial value *is* the value -- and
  # where each chain estimating its own showed up as chains that behaved
  # differently for no reason.
  fixed <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1,
    control = list(sampler = "nuts", chains = 3, warmup = 0, draws = 20, processes = FALSE,
      stepsize = 0.05))))
  expect_equal(fixed$sample$stepsize, rep(0.05, 3L))

  # Left unset, the chains estimate their own and need not agree; what is
  # asserted is only that the setting is not silently ignored.
  free <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1,
    control = list(sampler = "nuts", chains = 3, warmup = 0, draws = 20, processes = FALSE))))
  expect_false(isTRUE(all.equal(free$sample$stepsize, rep(0.05, 3L))))

  # And with warmup it is a starting point rather than the answer.
  moved <- suppressWarnings(suppressMessages(ctFitUncertainty(fit,
    uncertainty = "sample", cores = 1,
    control = list(sampler = "nuts", chains = 2, warmup = 60, draws = 20, processes = FALSE,
      stepsize = 0.05))))
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
# The joint-target sample with its effect draws saved, which the four tests
# below read and none writes into. It was drawn four times, identically -- the
# same fixture, settings and random number stream -- and is drawn once now,
# as the fixture is.
.sample_joint_effects <- function() fit_cached("sample_joint_effects",
  suppressWarnings(suppressMessages(ctFitUncertainty(.sample_fixture(),
    uncertainty = "sample", cores = 1,
    control = list(chains = 1, warmup = 60, draws = 60, saveEffects = TRUE,
      intoverpop = "none")))))

test_that("generation at given effects is exact and per-subject", {
  skip_without_julia()
  # intoverpop='none': each of these reads fit$sample$effects, which only exists
  # under the joint target -- the marginal default has no random effect in
  # the sampled vector to save.
  fit <- .sample_joint_effects()
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
  # intoverpop='none': each of these reads fit$sample$effects, which only exists
  # under the joint target -- the marginal default has no random effect in
  # the sampled vector to save.
  fit <- .sample_joint_effects()
  # Silence is the failure mode: the draws are summarised by default, so
  # without a message a user asking for a posterior predictive would get one
  # conditioned on point estimates and no way to notice.
  stripped <- fit
  stripped$sample$effects <- NULL
  # effects = 'fitted': the default draws new subjects and needs no draws.
  expect_message(ctGenerateFromFit(stripped, nsamples = 3,
    fullposterior = TRUE, cores = 1, effects = "fitted"), "saveEffects=TRUE")
  # With them, nothing to report.
  expect_no_message(ctGenerateFromFit(fit, nsamples = 3,
    fullposterior = TRUE, cores = 1, effects = "fitted"))
})

test_that("a sampled fit reports subject parameters from its own effects", {
  skip_without_julia()
  fit <- .sample_joint_effects()
  posterior <- fit$estimate$rawposterior
  draws <- ctSubjectPars(fit, pointest = FALSE)
  expect_equal(dim(draws)[1], nrow(posterior))

  # Draw s of a subject is its draw-s parameters plus its draw-s deviation:
  # on the raw scale exactly, so through a monotone transform the ranks over
  # draws agree perfectly. Modes linearised around the estimate would not.
  index <- fit$sample$effectIndex
  shared <- intersect(as.character(unique(index$parameter)),
    intersect(colnames(posterior), dimnames(draws)$param))
  expect_gt(length(shared), 0L)
  # Checked as monotonicity over the draws rather than a rank correlation of
  # exactly one: the engine rebuilds a draw's raw value through the population
  # factor, equal to `posterior + effects` only to rounding, so two draws whose
  # raw values agree to rounding may sort either way -- one such pair in 400
  # draws is a Spearman correlation of 1 - 2e-7. Pairs further apart than
  # rounding must all move the same way.
  for (p in shared) for (i in c(1L, max(index$subject))) {
    raw <- posterior[, p] + fit$sample$effects[, which(index$parameter == p &
      index$subject == i)]
    o <- order(raw)
    apart <- diff(raw[o]) > 1e-9 * max(1, abs(raw))
    steps <- sign(diff(draws[o, i, p]))[apart]
    expect_gt(sum(apart), length(raw) / 2)
    expect_equal(length(unique(steps[steps != 0])), 1L,
      info = paste(p, "subject", i))
  }

  # The same pairing through the point path: a fit whose estimate and mean
  # effects are draw s reports draw s.
  for (s in c(1L, nrow(posterior))) {
    one <- fit
    one$estimate$raw <- as.numeric(posterior[s, ])
    one$sample$effect_mean <- as.numeric(fit$sample$effects[s, ])
    expect_equal(as.numeric(ctSubjectPars(one)[1, , , drop = FALSE]),
      as.numeric(draws[s, , , drop = FALSE]))
  }

  # Without the effect draws the draws fall back to linearised modes, and say so.
  stripped <- fit
  stripped$sample$effects <- NULL
  expect_message(ctSubjectPars(stripped, pointest = FALSE), "saveEffects=TRUE")
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
  # intoverpop='none': each of these reads fit$sample$effects, which only exists
  # under the joint target -- the marginal default has no random effect in
  # the sampled vector to save.
  fit <- .sample_joint_effects()
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
  # intoverpop='none': each of these reads fit$sample$effects, which only exists
  # under the joint target -- the marginal default has no random effect in
  # the sampled vector to save.
  fit <- .sample_joint_effects()
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

test_that("a sampled fit's Hessian is not reused as curvature at its mean", {
  skip_without_julia()
  fit <- .sample_fixture()
  sampled <- suppressWarnings(suppressMessages(
    ctFitUncertainty(fit, uncertainty = "sample", cores = 1,
      control = list(chains = 1, warmup = 40, draws = 40))))

  # The exact Hessian a sampled fit carries was taken at the point the sampler
  # was placed from, and `$estimate$raw` is the posterior mean -- a different
  # point. `evaluated_at` is what says which, and it is what stops
  # `.ctBackendHessian()` reusing the matrix for the mean.
  at <- sampled$uncertainty$evaluated_at
  expect_false(is.null(at))
  # The fit's estimate, where its Hessian was taken -- not SAEM's placement
  # point, which the chains started from.
  expect_equal(as.numeric(at), as.numeric(fit$estimate$raw))
  # The premise: the two points really are different, or this guards nothing.
  expect_false(isTRUE(all.equal(as.numeric(at),
    as.numeric(sampled$estimate$raw), tolerance = 1e-8)))

  stored <- sampled$uncertainty$hessian
  # Asked about the point it was evaluated at, the cache answers with it.
  expect_identical(.ctBackendHessian(sampled, as.numeric(at)), stored)
  # Asked about the posterior mean, it must not: a guard written against
  # `$estimate$raw` instead of `evaluated_at` passes here and returns curvature
  # from the Laplace point, which is the defect this exists to prevent.
  atmean <- .ctBackendHessian(sampled, as.numeric(sampled$estimate$raw))
  expect_false(isTRUE(all.equal(atmean, stored, tolerance = 0)))

  # And an optimised fit's Hessian describes its estimate -- evaluated there,
  # or within a hundredth of a standard error where the finish kept the
  # Hessian it took at the hand-over -- so the same guard reuses it there,
  # which is the saving the field exists to keep.
  expect_lte(.ctBackendHessianDistance(fit$uncertainty$hessian,
    fit$uncertainty$evaluated_at, fit$estimate$raw), .ctBackendHessianReuse())
  expect_identical(.ctBackendHessian(fit, as.numeric(fit$estimate$raw)),
    fit$uncertainty$hessian)
})

# The placement pipeline: `ctFit(optimize = FALSE)` places the sampler through
# the same stages `optimize = TRUE` runs -- start, prior warm-up, substep
# mesh, the approach, the endgame's certification and its resume -- rather
# than a bare, uncertified optimisation (section 2f,
# review/OPTIM-consolidation-plan-2026-09-25.md). Before this, a start whose
# curvature had the wrong sign in some direction reached the sampler's metric
# exactly as it stood; certification is what a healthy `optimize = TRUE` fit
# always had and a sampled one never did.
test_that("a mixed-curvature start is certified before the sampler's metric is built", {
  skip_without_julia()
  # The fixture's data and model, not its fit: this test also needs the bare
  # `model_spec` (`fit = FALSE`) to verify its own precondition before trusting
  # what placement does with it.
  data <- laplace_fixture_data()
  model <- laplace_fixture_model()

  spec <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, intoverpop = "laplace", priors = TRUE,
    fit = FALSE)))
  npar <- .ctBackendNpar(spec)
  saddle <- rep(0, npar)
  saddle[1] <- 6

  # The precondition, checked directly rather than assumed: at this point the
  # exact Hessian genuinely has both signs present (a direction the objective
  # can still climb, alongside directions where it cannot), so certifying it
  # is not a no-op. `.ctBackendEndgameAt()` forms the Hessian at an arbitrary
  # point without optimising -- the same fit-free-assertion shape CLAUDE.md
  # asks for a matsetup column move, applied here to a starting point instead.
  before <- ctsem:::.ctBackendEndgameAt(spec, saddle)
  skip_if(is.null(before), "engine could not evaluate the endgame at this point")
  ev_before <- eigen(before$hessian, symmetric = TRUE, only.values = TRUE)$values
  skip_if_not(any(ev_before > 1e-6) && any(ev_before < -1e-6),
    "this start is no longer a mixed-curvature point on this engine version; pick a new one")

  # `placement = 'fit'`, the route this is about: under the default SAEM
  # placement the placement fit stops after the prior warm-up and SAEM moves
  # the start, so there is no certification to check.
  sampled <- suppressWarnings(suppressMessages(ctFit(data, model,
    backend = "julia", cores = 1, intoverpop = "laplace", priors = TRUE,
    optimize = FALSE, inits = saddle,
    sampleControl = list(chains = 1, warmup = 20, draws = 20, placement = "fit"))))

  # Certified, not merely run: the placement's own verdict, which is what
  # `$optim$converged` reports on a sampled fit now (decision 5,
  # review/OPTIM-consolidation-plan-2026-09-25.md) rather than the hardcoded
  # `TRUE` it used to be. `$optim$hessian_evaluations` is on the fit at all
  # only because the certification stage ran -- the bare `.ctJuliaOptimise()`
  # call this replaced had no such field to report.
  expect_true(isTRUE(sampled$optim$converged))
  expect_false(is.null(sampled$optim$hessian_evaluations))
  expect_gte(sampled$optim$hessian_evaluations, 1L)

  # The decisive check: the metric handed to the sampler has no direction the
  # objective can still climb. The starting point above did; this is what
  # "certified away before the metric is built" means, not merely that the
  # run finished.
  ev_after <- eigen(sampled$uncertainty$hessian, symmetric = TRUE,
    only.values = TRUE)$values
  expect_true(all(is.finite(ev_after)))
  expect_true(all(ev_after < 1e-6))
})
