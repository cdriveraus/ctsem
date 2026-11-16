# Random effects that reach dynamics depending on the state are linearised by
# the augmented route over their population spread, which can overstate the
# likelihood and move the estimate: on a milk-yield model whose yield drift
# reads fat, the augmented maximum read -2843 by its filter and about -3467 by
# a particle filter (scripts/nlgap-2026-11). ctFit() warns, on julia's
# intoverpop = 'augmented' and on stan, which has only that route. The check
# runs before the backend dispatch, so these use `fit = FALSE` and fit nothing.

.nlwarn_data <- function() {
  set.seed(1)
  data.frame(cow = rep(1:3, each = 5), day = rep(1:5, 3),
    yield = stats::rnorm(15), fat = stats::rnorm(15))
}

# Yield and fat; the drift of yield on fat reads fat. Every model below
# differs from it only in DRIFT and in which parameters vary; by default the
# initial and manifest means do.
.nlwarn_model <- function(DRIFT = matrix(c('drift_yield', 'cyf + nyf * fat',
    'drift_fat_yield', 'drift_fat'), 2, 2, byrow = TRUE), varying = NULL) {
  m <- suppressMessages(ctModel(type = 'ct', LAMBDA = diag(2),
    manifestNames = c('yield', 'fat'), latentNames = c('yield', 'fat'),
    PARS = matrix(c('cyf', 'nyf'), 2, 1), DRIFT = DRIFT,
    time = 'day', id = 'cow'))
  if (!is.null(varying)) m$pars$indvarying <- m$pars$param %in% varying &
      is.na(m$pars$value)
  m
}

# This warning alone: ctFit says other things about these models.
.nlwarn <- function(model, ...) {
  w <- character(0)
  withCallingHandlers(
    suppressMessages(ctFit(.nlwarn_data(), model, fit = FALSE, ...)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart('muffleWarning')
    })
  grep('dynamics that depend on the state', w, value = TRUE)
}

test_that('random initial means through a state dependent drift warn on stan', {
  w <- .nlwarn(.nlwarn_model(), backend = 'stan')
  expect_length(w, 1)
  expect_match(w, "T0m_yield, T0m_fat over")
  expect_match(w, "backend='julia'")
  # A random effect read by the state dependent cell itself.
  w <- .nlwarn(.nlwarn_model(varying = 'nyf'), backend = 'stan')
  expect_length(w, 1)
  expect_match(w, 'effects on nyf over')
})

test_that('linear dynamics and measurement effects alone are silent on stan', {
  linear <- matrix(c('drift_yield', 'cyf', 'drift_fat_yield', 'drift_fat'), 2, 2,
    byrow = TRUE)
  expect_length(.nlwarn(.nlwarn_model(linear), backend = 'stan'), 0)
  expect_length(.nlwarn(.nlwarn_model(varying = c('mm_yield', 'mm_fat')),
    backend = 'stan'), 0)
})

test_that('an effect the drift does not carry to the state read is silent', {
  # Fat's own drift reads fat, and the processes do not drive each other:
  # yield's effects cannot reach fat, fat's initial mean does. PARS[1,1] (cyf)
  # is read by yield's drift alone: it was once read into every transform
  # (`.ctContextCellTable()`), fat's initial mean among them.
  decoupled <- matrix(c('cyf', 0, 0, 'nyf * fat'), 2, 2, byrow = TRUE)
  expect_length(.nlwarn(.nlwarn_model(decoupled,
    varying = c('T0m_yield', 'cyf')), backend = 'stan'), 0)
  w <- .nlwarn(.nlwarn_model(decoupled, varying = c('T0m_yield', 'T0m_fat')),
    backend = 'stan')
  expect_length(w, 1)
  expect_match(w, 'effects on T0m_fat over')
})

test_that('julia warns on augmented only', {
  skip_without_julia()
  m <- .nlwarn_model()
  w <- .nlwarn(m, backend = 'julia', intoverpop = 'augmented')
  expect_length(w, 1)
  expect_match(w, "^intoverpop='augmented'")
  expect_match(w, 'The default intoverpop')
  expect_length(.nlwarn(m, backend = 'julia'), 0)
  expect_length(.nlwarn(m, backend = 'julia', optimize = FALSE, priors = TRUE), 0)
})
