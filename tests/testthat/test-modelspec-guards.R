# What a matrix cell means is decided by what its text looks like, so a cell
# that means one thing and reads as another builds a different model. These
# are the readings that used to happen silently, or that failed only in the
# generated program's own words. See R/ctModelSpecCheck.R.

.guard_model <- function(..., n.latent = 1, latentNames = 'e1') {
  args <- list(type = 'ct', n.latent = n.latent, n.manifest = 1,
    manifestNames = 'Y1', latentNames = latentNames, LAMBDA = matrix(1),
    Tpoints = 5)
  extra <- list(...)
  args[names(extra)] <- extra
  suppressMessages(do.call(ctModel, args))
}

test_that("a cell that is exactly a latent state name says so", {
  # Read as a name it is a free parameter; read as a reference it is the
  # state's value, and ctsem takes the second. On this model that is 11
  # parameters rather than 27, and nothing said which was built.
  expect_warning(.guard_model(MANIFESTMEANS = matrix('e1')),
    'latent state name')
  expect_warning(.guard_model(CINT = matrix('e1')), 'state\\[1\\]')
  # The warning names the explicit spelling, and that spelling is silent.
  expect_silent(suppressMessages(
    .guard_model(MANIFESTMEANS = matrix('state[1]'))))
  # An expression is not ambiguous and says nothing.
  expect_silent(suppressMessages(.guard_model(MANIFESTMEANS = matrix('e1*2'))))
  # An ordinary name is untouched.
  expect_silent(suppressMessages(.guard_model(MANIFESTMEANS = matrix('mm'))))
})

test_that("a cell that is exactly a tdpred name says so", {
  expect_warning(.guard_model(n.TDpred = 1, TDpredNames = 'tdp1',
    MANIFESTMEANS = matrix('tdp1')), 'time dependent predictor name')
})

test_that("an expression referencing nothing computable is refused", {
  # A compound cell becomes a calculation in the generated program, and what
  # makes it one is a reference to a state, a tdpred or a PARS cell. Without
  # one, nothing registers the names: the text became a free parameter's label
  # and its transform reached the backend as "99999 + 99999*".
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)')),
    'references no latent state')
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)')),
    'dr\\|-exp\\(param\\)')
  expect_error(.guard_model(DRIFT = matrix('dr*2')), 'nothing for it to be')
  expect_error(.guard_model(DRIFT = matrix('dr+dr2')), 'references no latent')
  # `dt` is not available to a cell, and the cell that looked like it scaled
  # with the interval did not.
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)*dt')),
    'references no latent')
  # Arithmetic on numbers is not an expression anyone means either.
  expect_error(.guard_model(DRIFT = matrix('2/3')), 'write the number itself')
})

test_that("a fresh parameter name inside a computable expression is fine", {
  # This is the working idiom and the reason the check is about the reference
  # rather than the names: `lbystate` is declared by being written here.
  expect_silent(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), Tpoints = 5,
    LAMBDA = matrix(c('lbystate * eta2 + 1', 0, 0, 1), 2, 2))))
  expect_silent(suppressMessages(.guard_model(n.TDpred = 1,
    TDpredNames = 'TD1', DRIFT = matrix('-log1p(exp(dr11)) * (1+TD1)'))))
  # And with the name declared in PARS, so is the expression that was refused
  # above.
  expect_silent(suppressMessages(.guard_model(DRIFT = matrix('-exp(dr)'),
    PARS = matrix('dr'))))
  expect_silent(suppressMessages(.guard_model(DRIFT = matrix('0.5*PARS[1,1]'),
    PARS = matrix('dr'))))
  expect_silent(suppressMessages(.guard_model(DRIFT = matrix('PARS[1,1]*1e-5'),
    PARS = matrix('dr'))))
  expect_silent(suppressMessages(.guard_model(n.TDpred = 1,
    TDpredNames = 'tdp1', MANIFESTMEANS = matrix('tdpreds[rowi, 1]'))))
})

test_that("an index outside what it indexes names the bound", {
  # Both of these reported "subscript out of bounds" from inside the matrix
  # unfolding, naming neither the cell nor the bound.
  expect_error(.guard_model(MANIFESTMEANS = matrix('state[9]')),
    'has 1 latent state')
  expect_error(.guard_model(DRIFT = matrix('PARS[2,2]'), PARS = matrix('dr')),
    'PARS is 1x1')
  expect_error(.guard_model(DRIFT = matrix('PARS[1,1]')), 'no PARS matrix')
})

test_that("sdscale that cannot apply says so", {
  # sdscale multiplies the prior on a population sd, so it means nothing for a
  # parameter with no population distribution. It was accepted and dropped.
  expect_warning(.guard_model(DRIFT = matrix('dr|||0.5')), 'does not vary')
  expect_warning(.guard_model(DRIFT = matrix('dr, sdscale=0.5')),
    'does not vary')
  expect_silent(suppressMessages(
    .guard_model(DRIFT = matrix('dr, indvarying=TRUE, sdscale=0.5'))))
  # MANIFESTMEANS varies by default, so the same field is meaningful there.
  expect_silent(suppressMessages(
    .guard_model(MANIFESTMEANS = matrix('mm, sdscale=0.5'))))
})

test_that("two cells asking for different sdscales for one parameter warn", {
  # Repeating a name is the equality constraint, so the cells are one
  # parameter with one prior scale. Which one was used depended on cell order.
  expect_warning(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    MANIFESTMEANS = matrix(c('mm||TRUE|0.5', 'mm||TRUE|3'), 2, 1))),
    'one prior scale')
})

test_that("an NA cell names the matrix and the cell", {
  # Four `[1]=="auto"` tests compare against the first cell, so a leading NA
  # came out as "missing value where TRUE/FALSE needed".
  expect_error(.guard_model(MANIFESTMEANS = matrix(NA)),
    'MANIFESTMEANS\\[1,1\\] is NA')
  expect_error(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    DRIFT = matrix(c('dr1', NA, '0', 'dr2'), 2, 2))), 'DRIFT\\[2,1\\] is NA')
})

test_that("a compound cell cannot carry fields, in either spelling", {
  # The message used to talk about `|` separators for a cell written with none.
  expect_error(.guard_model(MANIFESTMEANS = matrix('e1, indvarying=TRUE')),
    'cannot also carry transform, indvarying')
  expect_error(.guard_model(MANIFESTMEANS = matrix('e1||TRUE')),
    'cannot also carry transform, indvarying')
  expect_false(grepl('\\| separators',
    tryCatch(.guard_model(MANIFESTMEANS = matrix('e1, indvarying=TRUE')),
      error = function(e) conditionMessage(e))))
})

test_that("an empty tipreds field means no predictor acts on the cell", {
  # 'mm||||' has always meant that; the named form silently meant the default
  # instead, so the two spellings of one cell disagreed.
  expect_equal(.ctCellSpecToPipe('mm, tipreds='), 'mm||||')
  expect_equal(.ctCellSpecToPipe('mm, tipreds=c()'), 'mm||||')
  expect_equal(ctParSpec('mm', tipreds = character(0)), 'mm||||')
  # NULL is "not stated", which is a different cell: it keeps the default.
  expect_equal(ctParSpec('mm', tipreds = NULL), 'mm')

  eff <- function(spec) {
    m <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
      manifestNames = 'Y1', latentNames = 'e1', LAMBDA = matrix(1),
      n.TIpred = 1, TIpredNames = 'age', MANIFESTMEANS = matrix(spec)))
    m$pars$age_effect[m$pars$matrix %in% 'MANIFESTMEANS']
  }
  expect_equal(eff('mm||||'), 'FALSE')
  expect_equal(eff('mm, tipreds='), 'FALSE')
  expect_equal(eff('mm'), 'TRUE')

  # And the same through the assignment path, which reads the fields rather
  # than reparsing the cell. paste0() recycling made an empty list of
  # predictors look like a request for one with no name.
  m <- suppressMessages(ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
    manifestNames = 'Y1', latentNames = 'e1', LAMBDA = matrix(1),
    n.TIpred = 1, TIpredNames = 'age', MANIFESTMEANS = matrix('mm')))
  m$matrices$MANIFESTMEANS[1, 1] <- 'mm, tipreds='
  expect_equal(m$pars$age_effect[m$pars$matrix %in% 'MANIFESTMEANS'], 'FALSE')
  m$matrices$MANIFESTMEANS[1, 1] <- 'mm, tipreds=age'
  expect_equal(m$pars$age_effect[m$pars$matrix %in% 'MANIFESTMEANS'], 'TRUE')
})

test_that("the fields reader sees an empty trailing tipreds field", {
  # strsplit drops trailing empties, so the fifth field of 'mm||||' was
  # invisible and indistinguishable from not stating one.
  expect_equal(.ctCellSpecFields('mm||||')$tipreds, character(0))
  expect_null(.ctCellSpecFields('mm')$tipreds)
  expect_equal(.ctCellSpecFields('mm||||age')$tipreds, 'age')
})

test_that("a bare name for something a cell cannot reference says so", {
  # Unlike a latent or tdpred name, these stay ordinary free parameters --
  # which is the problem: the model builds and fits, and reports a parameter
  # called `age` beside a time independent predictor called `age`.
  expect_warning(.guard_model(n.TIpred = 1, TIpredNames = 'age',
    MANIFESTMEANS = matrix('age')), 'also a time independent predictor')
  expect_warning(.guard_model(MANIFESTMEANS = matrix('Y1')),
    'also a manifest variable')
  expect_warning(.guard_model(DRIFT = matrix('dt')),
    'time interval is not available')
  # A name that only collides with an identifier in the generated program is
  # harmless -- it is a label, and never emitted as code.
  expect_silent(suppressMessages(.guard_model(MANIFESTMEANS = matrix('state'))))
  expect_silent(suppressMessages(.guard_model(MANIFESTMEANS = matrix('PARS'))))
  expect_silent(suppressMessages(.guard_model(MANIFESTMEANS = matrix('DRIFT'))))
})

test_that("a refused expression names the unreachable thing it reached for", {
  # The error is already fatal, so saying which confusion this was costs
  # nothing and is the whole reason the word was written.
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)*dt')),
    'the time interval, which a cell cannot see')
  expect_error(.guard_model(n.TIpred = 1, TIpredNames = 'age',
    DRIFT = matrix('dr*age')), 'acts through the cell\'s tipreds field')
  expect_error(.guard_model(DRIFT = matrix('dr*Y1')),
    'a manifest variable, which a cell cannot see')
})
