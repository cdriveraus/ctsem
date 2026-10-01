# One effective sample size for every route that produces draws, and a warning
# from each when it ends short of it. Fit-free: the defaults, the sampler's
# budget rule and the three shortfall rules, on hand-built inputs.

test_that("every draw-producing route defaults to the same effective sample size", {
  target <- ctsem:::.ctEssTarget
  expect_identical(target, 200)
  # The exported corrections write the number out in their signatures.
  expect_identical(formals(ctsem::ctLaplaceCorrect)$target_ess, target)
  expect_identical(formals(ctsem::ctParticleCorrect)$target_ess, target)
  # The internal ones name it.
  expect_identical(eval(formals(ctsem:::.ctOptimImisDraws)$target_ess,
    asNamespace("ctsem")), target)
  expect_identical(eval(formals(ctsem:::imis_is)$target_ess,
    asNamespace("ctsem")), target)
  expect_identical(ctsem:::.ctBackendSampleControl(list())$min_ess, target)
  # uncertainty = 'is' sets isESS inside the draw stage when it is not given.
  body_text <- paste(deparse(body(ctsem:::.ctOptimDrawSamples)), collapse = "\n")
  expect_match(body_text, "control$isESS <- .ctEssTarget", fixed = TRUE)
})

test_that("the sampler's budget reaches past the draws asked for only with a target", {
  budget <- ctsem:::.ctBackendSampleBudget
  # A target: the first batch is sized from it, the budget is four times the
  # count asked for.
  b <- budget(500L, 4L, list(min_ess = 200))
  expect_identical(b$first, 50L)
  expect_identical(b$max_draws, 2000L)
  # Small targets keep the 50-draw floor; a count below it is respected.
  expect_identical(budget(500L, 1L, list(min_ess = 200))$first, 200L)
  expect_identical(budget(30L, 4L, list(min_ess = 200))$first, 30L)
  # maxDraws given is the budget, and the run starts with the count asked for.
  given <- budget(500L, 4L, list(min_ess = 200, max_draws = 900L))
  expect_identical(given$first, 500L)
  expect_identical(given$max_draws, 900L)
  # No target (minESS = 0): exactly the count asked for.
  none <- budget(500L, 4L, ctsem:::.ctBackendSampleControl(list(minESS = 0)))
  expect_identical(none$first, 500L)
  expect_null(none$max_draws)
})

test_that("a sample that ends short of its target says so, and one without a target keeps the old floor", {
  diag <- function(ess, target) list(chains = 4L, warmup = 200L, draws = 500L,
    divergent = 0L, rhat = c(a = 1.001, b = 1.002), ess = c(a = 900, b = ess),
    saturated = 0L, max_depth = 10L, unidentified = character(0),
    ess_target = target)
  floor <- ctsem:::.ctSampleEssFloor
  expect_identical(floor(diag(150, 200)), 200)
  expect_identical(floor(diag(150, NA_real_)), 100)
  expect_identical(floor(list()), 100)
  # Short of 200: warned, naming the target and the budget to raise.
  expect_warning(ctsem:::.ctSampleWarn(diag(150, 200)),
    "short of the target of 200.*maxDraws")
  expect_match(ctsem:::.ctSampleDiagnosis(diag(150, 200))$problems,
    "smallest effective sample size 150 \\(target 200\\)")
  # Met: nothing said, nothing recorded.
  expect_no_warning(ctsem:::.ctSampleWarn(diag(250, 200)))
  expect_length(ctsem:::.ctSampleDiagnosis(diag(250, 200))$problems, 0L)
  # No target: 150 passes, 80 does not, with the old advice.
  expect_no_warning(ctsem:::.ctSampleWarn(diag(150, NA_real_)))
  expect_warning(ctsem:::.ctSampleWarn(diag(80, NA_real_)),
    "set sampleControl\\$minESS")
  # Chains that disagree: with a target the run chose its own length, so the
  # advice is the budget, not the count asked for.
  rough <- function(target) utils::modifyList(diag(900, target),
    list(rhat = c(a = 1.001, b = 1.05)))
  expect_warning(ctsem:::.ctSampleWarn(rough(200)),
    "Raise sampleControl\\$maxDraws.*past the 500 draws per chain")
  expect_warning(ctsem:::.ctSampleWarn(rough(NA_real_)),
    "Raise the draw count -- iter in ctFit \\(now 700")
})

test_that("importance sampling that ends short of its target says so, not only below half of it", {
  report <- function(ess) ctsem:::.ctOptimImisReport(list(ess = ess), 200,
    weighted = TRUE, remedy = "REMEDY")
  # 150 against 200 passed silently when the bar was half the target.
  expect_warning(report(150), "150, short of its target of 200.*REMEDY")
  expect_no_warning(report(200))
  # The fixed-batch rule elsewhere has no target and does not claim one.
  expect_warning(ctsem:::.ctOptimEffectiveSampleWarn(40, floor = 50,
    remedy = "R", ndraws = 400), "40 from 400 draws\\. The intervals")
})
