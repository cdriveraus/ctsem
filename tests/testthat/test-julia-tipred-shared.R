# Constraining two parameters to one predictor effect, end to end.
#
# `<TI>_effect` may carry a *name*, and two parameters given the same name are
# meant to share one coefficient. For a while the name was parsed, stored and
# never read: `.ctTipredEffectLabel()` had no callers, and both backends
# allocated a fresh coefficient per parameter regardless, so a specification
# that asked for a constraint silently got none.
#
# The two backends still number coefficients independently -- stan as it walks
# the cells (R/ctModelWriter.R), julia predictor-major
# (`.ctJuliaTIEffects()`) -- so the indices differ by construction and what has
# to agree is which effects are *one* coefficient. That is what this file
# checks, on the spec and then on a fitted model.

.tishare_data <- function(n = 30) {
  set.seed(1)
  age <- stats::rnorm(n)
  do.call(rbind, lapply(seq_len(n), function(i) data.frame(
    id = i, time = 1:4,
    Y1 = stats::rnorm(4) + 2 * age[i], Y2 = stats::rnorm(4) + 2 * age[i],
    age = age[i])))
}

# Two manifest means, one predictor, everything else fixed: the smallest model
# in which two parameters can share an effect.
.tishare_model <- function(cells) {
  suppressWarnings(suppressMessages(ctModel(type = 'ct', n.latent = 1,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'), latentNames = 'e1',
    LAMBDA = matrix(c(1, 1), 2, 1), n.TIpred = 1, TIpredNames = 'age',
    tipredDefault = FALSE, DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTVAR = diag(0.5, 2), MANIFESTMEANS = matrix(cells, 2, 1))))
}

test_that('a shared effect name is one coefficient in the prepared spec', {
  skip_on_cran()
  d <- .tishare_data()

  shared <- suppressWarnings(suppressMessages(ctFit(dat = d,
    model = .tishare_model(c('a||||age=sh', 'b||||age=sh')),
    fit = FALSE, backend = 'julia')))
  separate <- suppressWarnings(suppressMessages(ctFit(dat = d,
    model = .tishare_model(c('a||||age', 'b||||age')),
    fit = FALSE, backend = 'julia')))

  sh <- as.data.frame(shared$ti_effects)
  sp <- as.data.frame(separate$ti_effects)
  # Both parameters carry an effect either way.
  expect_equal(sort(sh$parameter), sort(sp$parameter))
  expect_equal(nrow(sh), 2L)
  # Shared: one coefficient for the two of them. Separate: two.
  expect_equal(length(unique(sh$coefficient)), 1L)
  expect_equal(length(unique(sp$coefficient)), 2L)
})

test_that('the same name under two predictors is not one coefficient', {
  skip_on_cran()
  # A coefficient multiplies one predictor's values, so sharing is scoped per
  # predictor: constraining an age effect equal to a sex effect is a different
  # claim, and not what a repeated label says.
  d <- .tishare_data()
  # A TI predictor has to be constant within subject.
  d$sex <- as.numeric(d$id %% 2)
  model <- suppressWarnings(suppressMessages(ctModel(type = 'ct', n.latent = 1,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'), latentNames = 'e1',
    LAMBDA = matrix(c(1, 1), 2, 1), n.TIpred = 2,
    TIpredNames = c('age', 'sex'), tipredDefault = FALSE,
    DRIFT = matrix(-0.5), DIFFUSION = matrix(1), T0VAR = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTVAR = diag(0.5, 2),
    MANIFESTMEANS = matrix(c('a||||age=sh', 'b||||sex=sh'), 2, 1))))
  spec <- suppressWarnings(suppressMessages(ctFit(dat = d, model = model,
    fit = FALSE, backend = 'julia')))
  effects <- as.data.frame(spec$ti_effects)
  expect_equal(nrow(effects), 2L)
  expect_equal(length(unique(effects$coefficient)), 2L)
})

test_that('a shared effect costs one parameter and reports one value', {
  skip_on_cran()
  # The constraint has to reach the estimates, not just the spec: one fewer
  # free parameter, and the two parameters' reported effects identical rather
  # than merely close. Identical standard errors are what distinguishes a
  # shared coefficient from two that happen to agree.
  d <- .tishare_data()

  shared <- suppressWarnings(suppressMessages(ctFit(dat = d,
    model = .tishare_model(c('a||||age=sh', 'b||||age=sh')),
    backend = 'julia', cores = 1, verbose = 0)))
  separate <- suppressWarnings(suppressMessages(ctFit(dat = d,
    model = .tishare_model(c('a||||age', 'b||||age')),
    backend = 'julia', cores = 1, verbose = 0)))

  ssh <- summary(shared)
  ssp <- summary(separate)
  expect_equal(as.integer(ssp$npars) - as.integer(ssh$npars), 1L)

  # Reported per parameter, which is what makes the constraint visible: two
  # rows that agree exactly.
  tsh <- as.data.frame(ssh$tipreds)
  expect_equal(nrow(tsh), 2L)
  expect_equal(tsh$mean[1], tsh$mean[2])
  expect_equal(tsh$sd[1], tsh$sd[2])
  # Unconstrained, the same two rows differ.
  tsp <- as.data.frame(ssp$tipreds)
  expect_false(isTRUE(all.equal(tsp$mean[1], tsp$mean[2])))
})

test_that('a fixed effect claims no coefficient', {
  skip_on_cran()
  # An effect fixed to a value carries the value and is not estimated, so it
  # must not take a coefficient slot. `.ctJuliaTIEffects()` tested for an
  # *active* effect rather than a free one, which was inert only because
  # nothing reached it with a fixed effect: the fitting path refuses one
  # outright, and on the generation path every parnumber is NA. Both are other
  # functions' behaviour, so this pins the local test instead.
  spec <- .ctTipredEffectSpec(c('TRUE', '4.3', 'FALSE', 'myeffect'))
  expect_equal(.ctTipredEffectFree(spec), c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(.ctTipredEffectActive(spec), c(TRUE, TRUE, FALSE, TRUE))

  # Fitting a fixed effect is refused by name, which is what keeps the two
  # backends from disagreeing about whether it is estimated.
  fixedmodel <- suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 1, n.manifest = 1, manifestNames = 'Y1', latentNames = 'e1',
    LAMBDA = matrix(1), n.TIpred = 1, TIpredNames = 'age',
    tipredDefault = FALSE, DRIFT = matrix(-0.5), DIFFUSION = matrix(1),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTVAR = matrix(0.5),
    MANIFESTMEANS = matrix('mm||||age=0.7'))))
  d <- .tishare_data()
  expect_error(suppressWarnings(suppressMessages(ctFit(dat = d,
    model = fixedmodel, fit = FALSE, backend = 'julia'))),
    'fixed to a value')

  # Generating with one is what a fixed effect is for, and still works.
  set.seed(1)
  generated <- suppressWarnings(suppressMessages(ctGenerate(fixedmodel,
    n.subjects = 5, Tpoints = 4, backend = 'julia')))
  expect_equal(nrow(generated), 20L)
  expect_true('age' %in% colnames(generated))
})

test_that('ctFit names the fix when given an unconverted omx model', {
  # A ctModel(type='omx') object is a list of matrices, so the first model
  # field ctFit reads is a NULL inside an `&&`: the failure was "invalid
  # argument type", naming neither the argument nor the conversion.
  omx <- suppressMessages(ctModel(type = 'omx', Tpoints = 5, n.latent = 1,
    n.manifest = 1, manifestNames = 'Y1', latentNames = 'e1',
    LAMBDA = matrix(1)))
  err <- tryCatch(ctFit(datalong = .tishare_data(), model = omx, fit = FALSE),
    error = function(e) conditionMessage(e))
  expect_match(err, "ctModel\\(type='omx'\\)")
  expect_match(err, 'ctModelConvertOMX\\(model\\)')
  expect_false(grepl('invalid argument type', err, fixed = TRUE))
})
