# ctModelLavaan(), ctModelMplus(), ctModelFormula() and ctModelSyntax():
# fit-free. Each reader must build the model ctModel() builds from the same
# cells, and the writer's code must build the model it was given.

quiet <- function(expr) suppressMessages(expr)

# A predictor effect on a fixed cell means nothing, and tipredDefault = TRUE
# flags one there.
effects_of_parameters <- function(pars) {
  for (cl in grep('_effect$', names(pars), value = TRUE)) {
    pars[[cl]][is.na(pars$param)] <- 'FALSE'
  }
  pars
}

expect_same_model <- function(a, b) {
  expect_equal(effects_of_parameters(a$pars), effects_of_parameters(b$pars),
    ignore_attr = TRUE)
  expect_identical(a$latentNames, b$latentNames)
  expect_identical(a$manifestNames, b$manifestNames)
  expect_identical(a$TDpredNames, b$TDpredNames)
  expect_identical(a$TIpredNames, b$TIpredNames)
  expect_identical(a$manifesttype, b$manifesttype)
  expect_identical(a$fixedpopmeans, b$fixedpopmeans)
}

cell <- function(model, matrix, row, col = 1) {
  model$pars[model$pars$matrix == matrix & model$pars$row == row &
      model$pars$col == col, ]
}

# The code is run in an environment of its own: with settings to assign it
# names its model `model`.
roundtrip <- function(model) {
  for (style in c('lavaan', 'mplus', 'formula')) {
    code <- NULL
    expect_no_warning(code <- ctModelSyntax(model, style))
    back <- quiet(eval(parse(text = code), new.env()))
    expect_same_model(model, back)
    for (field in .ctSynSettings(model)) expect_identical(back[[field]],
      model[[field]])
  }
}

# The covariance a model's cells give, with its free cells at `free`.
implied <- function(model, matrix, free) {
  p <- model$pars[model$pars$matrix == matrix, ]
  cells <- matrix(0, max(p$row), max(p$col))
  cells[cbind(p$row, p$col)] <- ifelse(is.na(p$value), free, p$value)
  .ctCovForward(cells, model$covmattransform)
}

test_that('the three styles and ctModel() write the same default model', {
  m1 <- ctModel(LAMBDA = diag(2), manifestNames = c('y1', 'y2'),
    latentNames = c('eta1', 'eta2'), silent = TRUE)
  m2 <- ctModelLavaan('
    eta1 =~ y1   # a comment
    eta2 =~ y2
    eta1 ~ eta1 + eta2
    eta2 ~ eta1 +
      eta2')
  m3 <- ctModelMplus('
    eta1 BY y1;  ! a comment
    eta2 BY y2;
    eta1 eta2 ON eta1 eta2;')
  m4 <- ctModelFormula(y1 ~ eta1, y2 ~ eta2,
    d(eta1) ~ eta1 + eta2, d(eta2) ~ eta1 + eta2)
  expect_same_model(m1, m2)
  expect_same_model(m1, m3)
  expect_same_model(m1, m4)
  roundtrip(m1)
})

test_that('paths that are not written are zero, the rest follows ctModel()', {
  m <- ctModelLavaan('
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4
    eta1 ~ stress')
  # loadings: first fixed to 1, the others free
  expect_identical(cell(m, 'LAMBDA', 1, 1)$value, 1)
  expect_identical(cell(m, 'LAMBDA', 2, 1)$param, 'lambda_y2_eta1')
  expect_identical(cell(m, 'LAMBDA', 4, 1)$value, 0)
  # auto effects free, cross effects zero
  expect_identical(cell(m, 'DRIFT', 1, 1)$param, 'drift_eta1')
  expect_identical(cell(m, 'DRIFT', 1, 2)$value, 0)
  expect_identical(cell(m, 'DRIFT', 2, 1)$value, 0)
  # a variable that is no process is a time dependent predictor
  expect_identical(m$TDpredNames, 'stress')
  expect_identical(cell(m, 'TDPREDEFFECT', 1, 1)$param, 'td_eta1_stress')
  expect_identical(cell(m, 'TDPREDEFFECT', 2, 1)$value, 0)
  # intercepts of indicators free and varying, of processes zero
  expect_true(cell(m, 'MANIFESTMEANS', 1)$indvarying)
  expect_identical(cell(m, 'CINT', 1)$value, 0)
  # noise covariance free, residual covariance zero
  expect_identical(cell(m, 'DIFFUSION', 2, 1)$param, 'diff_eta2_eta1')
  expect_identical(cell(m, 'MANIFESTVAR', 2, 1)$value, 0)
  roundtrip(m)
})

test_that('modifiers fix, free, label and constrain', {
  m <- ctModelLavaan('
    eta1 =~ NA*y1 + 0.8*y2 + l*y3
    eta2 =~ y4 + l*y5
    eta1 ~ a*eta1 + eta2
    eta2 ~ b*eta2
    a == b
    eta1 ~~ 0*eta2
    eta1 ~~ 0.25*eta1
    y1 ~~ y2
    y1 ~ 0*1
    T0m_eta1 ~ 0*1
    T0var_eta2_eta1 ~ 0*1')
  expect_identical(cell(m, 'LAMBDA', 1, 1)$param, 'lambda_y1_eta1')
  expect_identical(cell(m, 'LAMBDA', 2, 1)$value, 0.8)
  expect_identical(cell(m, 'LAMBDA', 3, 1)$param, 'l')
  expect_identical(cell(m, 'LAMBDA', 5, 2)$param, 'l')
  expect_identical(cell(m, 'DRIFT', 1, 1)$param, 'a')
  expect_identical(cell(m, 'DRIFT', 2, 2)$param, 'a')
  expect_identical(cell(m, 'DIFFUSION', 2, 1)$value, 0)
  # a variance of 0.25 is a standard deviation of 0.5
  expect_identical(cell(m, 'DIFFUSION', 1, 1)$value, 0.5)
  expect_identical(cell(m, 'MANIFESTVAR', 2, 1)$param, 'mvary2_y1')
  expect_identical(cell(m, 'MANIFESTMEANS', 1)$value, 0)
  expect_identical(cell(m, 'T0MEANS', 1)$value, 0)
  expect_identical(cell(m, 'T0VAR', 2, 1)$value, 0)
  roundtrip(m)

  mp <- ctModelMplus('
    eta1 BY y1* y2@0.8
            y3 (l);
    eta2 BY y4
            y5 (l);
    eta1 ON eta1 (a)
            eta2;
    eta2 ON eta2 (a);
    eta1 WITH eta2@0;
    eta1@0.25;
    y1 WITH y2;
    [y1@0];
    [T0m_eta1@0];
    [T0var_eta2_eta1@0];')
  expect_same_model(m, mp)

  mf <- ctModelFormula(
    y1 ~ 0 + NA*eta1, y2 ~ 0.8*eta1, y3 ~ "l"*eta1,
    y4 ~ eta2, y5 ~ "l"*eta2,
    d(eta1) ~ "a"*eta1 + eta2, d(eta2) ~ "a"*eta2,
    diff_eta2_eta1 ~ 0, diff_eta1 ~ 0.5, mvary2_y1 ~ 1,
    T0m_eta1 ~ 0, T0var_eta2_eta1 ~ 0)
  expect_same_model(m, mf)

  # in a formula 1 estimates, so a value of 1 is written 1*1
  one <- ctModelLavaan('
    eta1 =~ y1
    eta1 ~ 1*1
    y1 ~ 1*1
    T0m_eta1 ~ 1*1')
  expect_identical(cell(one, 'CINT', 1)$value, 1)
  expect_identical(cell(one, 'MANIFESTMEANS', 1)$value, 1)
  expect_identical(cell(one, 'T0MEANS', 1)$value, 1)
  roundtrip(one)
})

test_that('random effects and predictors of parameters', {
  viaModel <- quiet(ctModel(LAMBDA = diag(2), manifestNames = c('y1', 'y2'),
    latentNames = c('eta1', 'eta2'), TIpredNames = c('age', 'sex'),
    tipredDefault = FALSE, silent = TRUE,
    DRIFT = matrix(c('drift_eta1', 'a||TRUE||age,sex=b', 0, 'drift_eta2'), 2,
      byrow = TRUE),
    MANIFESTMEANS = c('mm_y1||FALSE', 'mm_y2||||age'),
    CINT = c('cint_eta1', 0)))
  ml <- ctModelLavaan('
    eta1 =~ y1
    eta2 =~ y2
    eta1 ~ eta1 + a*eta2
    eta1 ~ 1
    a ~ age + b*sex
    a ~~ a
    level: 2
    y1 ~~ 0*y1
    y2 ~ age')
  mp <- ctModelMplus('
    eta1 BY y1;
    eta2 BY y2;
    a | eta1 ON eta2;
    [eta1];
    %BETWEEN%
    a ON age
         sex (b);
    y1@0;
    y2 ON age;')
  mf <- ctModelFormula(
    y1 ~ eta1 + (0 | id),
    y2 ~ eta2,
    d(eta1) ~ 1 + eta1 + "a"*eta2 + (eta2 | id),
    a ~ age + "b"*sex,
    mm_y2 ~ age)
  expect_same_model(viaModel, ml)
  expect_same_model(viaModel, mp)
  expect_same_model(viaModel, mf)
  roundtrip(viaModel)
  # a parameter formula alone does the same as the random effect term
  mf2 <- ctModelFormula(y1 ~ eta1, y2 ~ eta2,
    d(eta1) ~ 1 + eta1 + eta2,
    drift_eta1_eta2 ~ "a"*1 + age + "b"*sex + (1 | id),
    mm_y1 ~ (0 | id), mm_y2 ~ age)
  expect_same_model(viaModel, mf2)
  # a fixed predictor effect
  m <- ctModelLavaan('
    eta1 =~ y1
    drift_eta1 ~ 0.3*age')
  expect_identical(cell(m, 'DRIFT', 1, 1)$age_effect, '0.3')
  roundtrip(m)
})

test_that('smooth effects of a predictor', {
  m <- ctModelFormula(y1 ~ eta1, d(eta1) ~ 1 + eta1,
    cint_eta1 ~ s(age, k = 6) + sex)
  expect_identical(names(m$TIsmooths), 'age')
  expect_identical(cell(m, 'CINT', 1)$age_lin_effect, 'TRUE')
  expect_identical(cell(m, 'CINT', 1)$sex_effect, 'TRUE')
  expect_identical(cell(m, 'DRIFT', 1, 1)$age_lin_effect, 'FALSE')
  ml <- ctModelLavaan('
    eta1 =~ y1
    eta1 ~ 1
    cint_eta1 ~ s(age, k = 6) + sex')
  mp <- ctModelMplus('
    eta1 BY y1;
    [eta1];
    cint_eta1 ON s(age, k = 6) sex;')
  expect_same_model(m, ml)
  expect_same_model(m, mp)
  roundtrip(m)
})

test_that('nested levels', {
  viaModel <- quiet(ctModel(LAMBDA = diag(1), manifestNames = 'y1',
    latentNames = 'eta1', id = c('person', 'school'), silent = TRUE,
    MANIFESTMEANS = 'mm_y1||school',
    DRIFT = 'drift_eta1||person'))
  ml <- ctModelLavaan(id = c('person', 'school'), '
    eta1 =~ y1
    level: person
    drift_eta1 ~~ drift_eta1
    level: school
    y1 ~~ y1')
  mp <- ctModelMplus(id = c('person', 'school'), '
    eta1 BY y1;
    %BETWEEN person%
    drift_eta1;
    %BETWEEN school%
    y1;')
  mf <- ctModelFormula(id = c('person', 'school'),
    y1 ~ eta1 + (1 | school),
    d(eta1) ~ eta1 + (eta1 | person))
  expect_same_model(viaModel, ml)
  expect_same_model(viaModel, mp)
  expect_same_model(viaModel, mf)
  roundtrip(viaModel)
  expect_error(ctModelLavaan(id = c('person', 'school'), '
    eta1 =~ y1
    level: class
    y1 ~~ y1'), 'not a level')
})

test_that('ctPar() and quoted cells give the other fields', {
  viaModel <- quiet(ctModel(LAMBDA = diag(2), manifestNames = c('y1', 'y2'),
    latentNames = c('eta1', 'eta2'), silent = TRUE, PARS = 'dr',
    DRIFT = matrix(c('drift_eta1', 'a|exp(param)|TRUE|0.5', '-exp(dr)*eta1',
      'drift_eta2'), 2, byrow = TRUE),
    CINT = c('lev=0||TRUE', 0)))
  ml <- ctModelLavaan(PARS = 'dr', '
    eta1 =~ y1
    eta2 =~ y2
    eta1 ~ eta1 + ctPar(a, transform = "exp(param)", indvarying = TRUE, sdscale = 0.5)*eta2
    eta2 ~ "-exp(dr)*eta1"*eta1 + eta2
    eta1 ~ ctPar("lev", mean = 0, indvarying = TRUE)*1')
  mp <- ctModelMplus(PARS = 'dr', '
    eta1 BY y1;
    eta2 BY y2;
    eta1 ON eta1
            eta2 (a|exp(param)|TRUE|0.5);
    eta2 ON eta1 ("-exp(dr)*eta1")
            eta2;
    [eta1] (lev, mean=0, indvarying=TRUE);')
  mf <- ctModelFormula(PARS = 'dr',
    y1 ~ eta1, y2 ~ eta2,
    d(eta1) ~ ctPar("lev", mean = 0, indvarying = TRUE)*1 + eta1 +
      ctPar("a", transform = "exp(param)", sdscale = 0.5)*eta2 + (eta2 | id),
    d(eta2) ~ "-exp(dr)*eta1"*eta1 + eta2)
  expect_same_model(viaModel, ml)
  expect_same_model(viaModel, mp)
  expect_same_model(viaModel, mf)
  roundtrip(viaModel)
  # random effects around a stated mean, the lme4 way
  m <- ctModelFormula(y1 ~ eta1, d(eta1) ~ eta1 + (1 | id))
  expect_identical(unname(m$fixedpopmeans['cint_eta1']), 0)
  expect_true(cell(m, 'CINT', 1)$indvarying)
  # and the Mplus way: the mean of a random slope fixed in %BETWEEN%
  s <- ctModelMplus('
    eta1 BY y1;
    eta2 BY y2;
    s | eta1 ON eta2;
    %BETWEEN%
    [s@0.2];
    s;')
  expect_identical(cell(s, 'DRIFT', 1, 2)$param, 's')
  expect_true(cell(s, 'DRIFT', 1, 2)$indvarying)
  expect_identical(unname(s$fixedpopmeans['s']), 0.2)
  roundtrip(s)
  # statements may come in any order
  expect_same_model(
    ctModelLavaan('a ~ age\n a ~~ a\n eta1 =~ y1\n eta2 =~ y2\n eta1 ~ a*eta2'),
    ctModelLavaan('eta1 =~ y1\n eta2 =~ y2\n eta1 ~ a*eta2\n a ~~ a\n a ~ age'))
})

test_that('indicator types', {
  m <- ctModelLavaan('
    eta1 =~ y1 + y2 + y3 + y4
    y2 | t1
    y3 | t1 + t2 + t3',
    manifesttype = c(y4 = 'count'))
  expect_identical(as.integer(m$manifesttype), c(0L, 1L, 2L, 3L))
  expect_identical(as.integer(m$ncategories[3]), 4L)
  mp <- ctModelMplus('
    eta1 BY y1-y4;
    [y2$1];
    [y3$1 y3$2 y3$3];',
    manifesttype = c(y4 = 'count'))
  expect_same_model(m, mp)
  mf <- ctModelFormula(y1 + y2 + y3 + y4 ~ eta1,
    family = list(y2 = 'bernoulli', y3 = 'cumulative', y4 = stats::poisson),
    ncategories = c(y3 = 4))
  expect_same_model(m, mf)
  roundtrip(m)
  mc <- ctModelFormula(y1 ~ eta1, manifesttype = 'censored',
    censormin = c(y1 = 0))
  expect_identical(as.integer(mc$manifesttype), 4L)
  expect_identical(mc$censormin, 0)
  roundtrip(mc)
})

test_that('discrete time, processes without indicators, variable order', {
  expect_message(m <- ctModelLavaan(type = 'dt', '
    eta1 =~ y1
    trend ~ 0*trend
    eta1 ~ eta1 + 1*trend'), 'trend is a process without indicators')
  expect_false(m$continuoustime)
  expect_identical(m$latentNames, c('eta1', 'trend'))
  expect_identical(cell(m, 'DRIFT', 2, 2)$value, 0)
  expect_identical(cell(m, 'DRIFT', 1, 2)$value, 1)
  roundtrip(m)
  # the writer keeps an order that is not the order of first appearance
  viaModel <- ctModel(LAMBDA = matrix(c(0, 1, 1, 0), 2), silent = TRUE,
    manifestNames = c('y1', 'y2'), latentNames = c('eta1', 'eta2'))
  roundtrip(viaModel)
})

test_that('a whole Mplus input, lags and ranges', {
  expect_message(m <- ctModelMplus('
    TITLE: an example;
    VARIABLE: NAMES = id time y1 y2 y3;
      CLUSTER = id;
    ANALYSIS: TYPE = TWOLEVEL RANDOM;
    MODEL:
      %WITHIN%
      f BY y1-y3 (&1);
      f ON f&1;
      %BETWEEN%
      y1-y3 ON age;
    OUTPUT: TECH1;', type = 'dt'), 'MODEL: section')
  expect_identical(m$manifestNames, c('y1', 'y2', 'y3'))
  expect_identical(cell(m, 'DRIFT', 1, 1)$param, 'drift_f')
  expect_identical(cell(m, 'MANIFESTMEANS', 3)$age_effect, 'TRUE')
  expect_error(ctModelMplus('f BY y1; f ON f&2;'), 'ctModelHigherOrder')
})

test_that('a lag is discrete time, and a value keeps its meaning', {
  # An autoregression of 0.8 is 0.8 over one occasion in the model it is read
  # into, and is not read as a continuous time auto effect of +0.8.
  m <- ctModelMplus('f BY y1; g BY y2; f ON f&1@0.8 g&1; g ON g&1;', type = 'dt')
  expect_false(m$continuoustime)
  expect_identical(cell(m, 'DRIFT', 1, 1)$value, 0.8)
  expect_identical(cell(m, 'DRIFT', 1, 2)$param, 'drift_f_g')
  expect_error(ctModelMplus('f BY y1; f ON f&1@0.8;'), "type = 'dt'")
  expect_error(ctModelMplus('f BY y1; f ON f&1;'), "type = 'dt'")
  expect_error(ctModelMplus('f BY y1; f&1 WITH f;', type = 'dt'), 'after ON')
  # the writer writes the lag where the model has one
  expect_match(ctModelSyntax(m, 'mplus'), 'f ON f&1@0.8 g&1;', fixed = TRUE)
  expect_false(grepl('&1', ctModelSyntax(ctModelMplus('f BY y1; f ON f;'),
    'mplus'), fixed = TRUE))
  roundtrip(m)
})

test_that('a covariance that is zero in the syntax is zero in the model', {
  # separate blocks: zero whatever the free cells are
  m <- ctModelLavaan('
    a =~ y1
    b =~ y2
    c =~ y3
    d =~ y4
    a ~~ 0*c + 0*d
    b ~~ 0*c + 0*d
    y1 ~~ y2')
  for (free in list(0.7, c(1.3, -0.4, 2))) {
    S <- implied(m, 'DIFFUSION', free)
    expect_equal(S[3:4, 1:2], matrix(0, 2, 2))
    expect_true(all(abs(S[2, 1]) > 0.01, abs(S[4, 3]) > 0.01))
    R <- implied(m, 'MANIFESTVAR', free)
    expect_equal(R[3:4, 1:2], matrix(0, 2, 2))
    expect_equal(R[4, 3], 0)
    expect_gt(abs(R[2, 1]), 0.01)
  }
  roundtrip(m)
  # a zero beside covariances that connect its two variables is refused,
  # written out or left out, in each matrix and style
  expect_error(ctModelLavaan('a =~ y1\n b =~ y2\n c =~ y3\n a ~~ 0*b'),
    'noise covariance of a and b')
  expect_error(ctModelLavaan('f =~ y1 + y2 + y3\n y1 ~~ y2\n y2 ~~ y3'),
    'residual covariance of y1 and y3')
  expect_error(ctModelMplus('a BY y1; b BY y2; c BY y3; a WITH b@0;'),
    'through c')
  expect_error(ctModelFormula(y1 ~ a, y2 ~ b, y3 ~ c, T0var_b_a ~ 0),
    'initial covariance of a and b')
  # what the refused model would have been: not a zero covariance
  cells <- diag(3)
  cells[3, 1:2] <- 0.5
  expect_gt(abs(.ctCovForward(cells, 'z')[2, 1]), 0.1)
  # the cell itself can still be set, in quotes, and a model that has such
  # cells is written that way
  q <- ctModelLavaan('a =~ y1\n b =~ y2\n c =~ y3\n a ~~ "0"*b')
  expect_identical(cell(q, 'DIFFUSION', 2, 1)$value, 0)
  expect_match(ctModelSyntax(q, 'mplus'), 'b WITH a ("0");', fixed = TRUE)
  roundtrip(q)
  viaModel <- ctModel(LAMBDA = diag(3), silent = TRUE,
    MANIFESTVAR = matrix(c('v1', 'v21', 0, 0, 'v2', 'v32', 0, 0, 'v3'), 3))
  roundtrip(viaModel)
})

test_that('a label names a first loading and NA frees it', {
  expect_message(m1 <- ctModelLavaan('f =~ a*y1 + y2'), 'stays fixed at 1')
  expect_message(m2 <- ctModelMplus('f BY y1 (a)\n y2;'), 'stays fixed at 1')
  expect_message(m3 <- ctModelFormula(y1 ~ "a"*f, y2 ~ f), 'stays fixed at 1')
  for (m in list(m1, m2, m3)) {
    expect_identical(cell(m, 'LAMBDA', 1, 1)$value, 1)
    expect_identical(cell(m, 'LAMBDA', 2, 1)$param, 'lambda_y2_f')
  }
  for (m in list(
    ctModelLavaan('f =~ NA*y1 + a*y1 + y2\n f ~~ 1*f'),
    ctModelMplus('f BY y1* (a)\n y2; f@1;'),
    ctModelFormula(y1 ~ NA*f + "a"*f, y2 ~ f, diff_f ~ 1*1))) {
    expect_identical(cell(m, 'LAMBDA', 1, 1)$param, 'a')
    expect_identical(cell(m, 'DIFFUSION', 1, 1)$value, 1)
    roundtrip(m)
  }
  # a label on any other loading labels it, free
  expect_identical(cell(ctModelLavaan('f =~ y1 + a*y2'), 'LAMBDA', 2, 1)$param,
    'a')
})

test_that('the writer carries the settings of a model, and says what it cannot', {
  m <- ctModelLavaan('
    eta1 =~ y1
    eta2 =~ y2
    eta2 ~ eta1 + eta2
    drift_eta2 ~ age')
  m$covmattransform <- 'cholesky'
  m$tipredeffectscale <- 3
  m$rawpopsdbase <- 'normal(0,0.5)'
  m$RAWPOPVAR[2, 1] <- '0'
  code <- ctModelSyntax(m, 'lavaan')
  expect_match(code, 'model$covmattransform <- "cholesky"', fixed = TRUE)
  expect_match(code, 'model$tipredeffectscale <- 3', fixed = TRUE)
  roundtrip(m)
  # a model at its defaults is still one call
  plain <- ctModelLavaan('eta1 =~ y1')
  expect_false(grepl('model', ctModelSyntax(plain, 'lavaan'), fixed = TRUE))
  # the comparison behind the warning names the cell and the field
  other <- plain
  other$pars$transform[other$pars$matrix == 'DRIFT'] <- 'exp(param)'
  other$tipredsimputedscale <- 2
  expect_identical(.ctSynDifferences(plain, plain), character())
  expect_identical(.ctSynDifferences(plain, other),
    c('transform of DRIFT[1,1]', 'tipredsimputedscale'))
})

test_that('a mistyped parameter name is not read as a new variable', {
  expect_error(ctModelFormula(y ~ f, drfit_f ~ age), 'close to the parameter name drift_f')
  expect_error(ctModelLavaan('f =~ y\n drfit_f ~ age'), 'drift_f')
  expect_error(ctModelMplus('f BY y; [T0M_f@0];'), 'T0m_f')
  # a variable that is named so is declared
  m <- ctModelFormula(y ~ f, drfit_f ~ f, manifestNames = c('y', 'drfit_f'))
  expect_identical(m$manifestNames, c('y', 'drfit_f'))
})

test_that('a model with default predictor effects is written out', {
  m <- quiet(ctModel(LAMBDA = diag(2), manifestNames = c('y1', 'y2'),
    latentNames = c('eta1', 'eta2'), TDpredNames = 'td', TIpredNames = 'age',
    silent = TRUE))
  roundtrip(m)
})

test_that('what the syntax cannot say is refused with a pointer', {
  expect_error(ctModelLavaan('eta1 =~ y1\n ab := a*b'), 'ctSummaryMatrices')
  expect_error(ctModelLavaan('eta1 =~ y1\n a > 0'), 'transform')
  expect_error(ctModelLavaan('group: 1\n eta1 =~ y1'), 'time independent')
  expect_error(ctModelLavaan('eta1 =~ y1\n eta1 ~ c(a, b)*eta1'), 'group')
  expect_error(ctModelLavaan('eta1 =~ y1\n y1 ~ age'), 'mm_y1')
  expect_error(ctModelLavaan('eta1 =~ y1 + y2\n y1 ~~ 0.3*y2'), 'ctCov')
  expect_error(ctModelLavaan('eta1 =~ y1\n nopar ~~ nopar'), 'covariance')
  expect_error(ctModelLavaan('eta1 =~ y1\n level: 2\n nopar ~ age'),
    'not a parameter')
  expect_error(ctModelLavaan('eta1 =~ y1\n y1 =~ y2'), 'process and as an')
  expect_error(ctModelLavaan('eta1 =~ y1\n eta1 ~ eta1*eta1'), 'variable')
  expect_error(ctModelLavaan('eta1 =~ y1', DRIFT = 'auto'), 'written by')
  expect_error(ctModelFormula(y1 ~ eta1, d(eta1) ~ eta1 + eta1:age),
    'interactions')
  expect_error(ctModelFormula(y1 ~ eta1, d(eta1) ~ a*eta1), 'quotes')
  expect_error(ctModelFormula(y1 ~ eta1, d(eta1) ~ eta1 + (1 | subject)),
    'not a level')
  expect_error(ctModelMplus('f BY y1; f XWITH f;'), 'XWITH')
  # a factor model between subjects, and uncorrelated random effects
  expect_error(ctModelLavaan('eta =~ y1 + y2\n level: 2\n fb =~ y1 + y2'),
    'popblocks')
  expect_error(ctModelMplus('eta BY y1 y2; %BETWEEN% fb BY y1 y2;'), 'popblocks')
  expect_error(ctModelFormula(y1 ~ eta, d(eta) ~ eta + (1 || id)), 'popunique')
  expect_message(ctModelLavaan('eta1 =~ y1\n eta1 ~ start(-0.5)*eta1'),
    'inits')
  # a mistyped field name in a quoted cell is ctPar's error
  expect_error(ctModelLavaan('eta1 =~ y1\n eta1 ~ "a, indvaring=TRUE"*eta1'),
    'indvarying')
})
