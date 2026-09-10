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

test_that("in PARS a bare state or tdpred name is refused, not warned about", {
  # A PARS cell exists to declare a parameter, so a bare name there reads as
  # the declaration. Taking it as a state reference instead is not a model
  # anyone wants, and the label would collide with the state name everywhere
  # else too.
  expect_error(.guard_model(PARS = matrix('e1'),
    DRIFT = matrix('-0.5 * (1 + PARS[1,1])')), 'latent state name')
  expect_error(.guard_model(n.TDpred = 1, TDpredNames = 'tdp1',
    PARS = matrix('tdp1'), DRIFT = matrix('-0.5 * (1 + PARS[1,1])')),
    'time dependent predictor name')
  # An expression over the same state is the supported way to write a
  # state-dependent PARS cell.
  expect_silent(suppressMessages(.guard_model(PARS = matrix('log1p(exp(e1))'),
    DRIFT = matrix('-0.5 * (1 + PARS[1,1])'))))
  # And an ordinary label is untouched.
  expect_silent(suppressMessages(.guard_model(PARS = matrix('p1'),
    DRIFT = matrix('-0.5 * (1 + PARS[1,1])'))))
})

test_that("an expression over undeclared names is refused", {
  # Nothing here was declared, so the cell became a free parameter whose label
  # was the expression text and whose transform reached the backend as
  # "99999 + 99999*".
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)')), 'not declared')
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)')),
    'dr\\|-exp\\(param\\)')
  expect_error(.guard_model(DRIFT = matrix('dr*2')), 'not declared')
  expect_error(.guard_model(DRIFT = matrix('dr+dr2')), '"dr", "dr2"')
  # Arithmetic on numbers has no names to report as undeclared, so it gets its
  # own message.
  expect_error(.guard_model(DRIFT = matrix('2/3')), 'names nothing at all')
  expect_error(.guard_model(DRIFT = matrix('2/3')), 'write the number itself')
})

test_that("a name used in an expression must be declared somewhere", {
  # Writing a name inside an expression does not declare it. Measured: with no
  # PARS cell for it, `-exp(zz)*eta1` dies at ctFit() with `object 'zz' not
  # found`, and so do the same shapes in LAMBDA, DIFFUSION, MANIFESTMEANS,
  # CINT and TDPREDEFFECT. The exception is worse rather than better --
  # `-log1p(exp(zz))*(1+TD1)` builds and `zz` is not a parameter at all -- so
  # none of these are models anyone can fit.
  expect_error(.guard_model(DRIFT = matrix('-exp(zz)*e1')), 'not declared')
  expect_error(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), Tpoints = 5,
    LAMBDA = matrix(c('lbystate * eta2 + 1', 0, 0, 1), 2, 2))),
    'not declared')
  expect_error(.guard_model(n.TDpred = 1, TDpredNames = 'TD1',
    DRIFT = matrix('-log1p(exp(dr11)) * (1+TD1)')), 'not declared')
  # The message says how, and names the parameter it would declare.
  msg <- tryCatch(.guard_model(DRIFT = matrix('-exp(zz)*e1')),
    error = function(e) conditionMessage(e))
  expect_match(msg, "PARS = c\\('zz'\\)")
})

test_that("declared in PARS, the same expressions are accepted", {
  expect_silent(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), Tpoints = 5,
    PARS = matrix('lbystate'),
    LAMBDA = matrix(c('lbystate * eta2 + 1', 0, 0, 1), 2, 2))))
  expect_silent(suppressMessages(.guard_model(n.TDpred = 1,
    TDpredNames = 'TD1', PARS = matrix('dr11'),
    DRIFT = matrix('-log1p(exp(dr11)) * (1+TD1)'))))
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
  # Some undeclared names are not names a PARS cell could supply, and saying
  # so is the point: "Add PARS = c('dt')" would send the reader the wrong way
  # entirely, because there is no parameter to declare.
  expect_error(.guard_model(DRIFT = matrix('-exp(dr)*dt')),
    'the time interval, which a cell cannot see')
  expect_error(.guard_model(n.TIpred = 1, TIpredNames = 'age',
    DRIFT = matrix('dr*age')), 'acts through the cell\'s tipreds field')
  expect_error(.guard_model(DRIFT = matrix('dr*Y1')),
    'a manifest variable, which a cell cannot see')
  # The PARS suggestion names only the parameter it could actually declare.
  msg <- tryCatch(.guard_model(DRIFT = matrix('-exp(dr)*dt')),
    error = function(e) conditionMessage(e))
  expect_match(msg, "PARS = c\\('dr'\\)")
  expect_false(grepl("PARS = c('dt')", msg, fixed = TRUE))
})

test_that("one parameter cannot be transformed two ways", {
  # Repeating a name is the equality constraint: one raw value, one random
  # effect. The transform stayed per cell, so that one value came out as two
  # different numbers and there was no answer to what the parameter was.
  mm2 <- function(a, b) suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 2, n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    MANIFESTMEANS = matrix(c(a, b), 2, 1))))

  # Two transforms written out.
  expect_error(mm2('mm|exp(param)', 'mm|log(param)'), 'transform it differently')
  # One written, one left at the cell default -- just as ambiguous.
  expect_error(mm2('mm|exp(param)', 'mm'), 'transform it differently')
  # The message names both cells and both transforms, and points at PARS.
  msg <- tryCatch(mm2('mm|exp(param)', 'mm|log(param)'),
    error = function(e) conditionMessage(e))
  expect_match(msg, 'MANIFESTMEANS\\[1,1\\] as exp\\(param\\)')
  expect_match(msg, 'MANIFESTMEANS\\[2,1\\] as log\\(param\\)')
  # The remedy names the parameter, which is what a user writes in a cell --
  # not the PARS coordinates, which ctsem substitutes for them -- and shows
  # the per-cell expression, which is what lets one parameter serve cells
  # needing different transforms.
  expect_match(msg, "PARS = c\\('mm'\\)")
  # The example carries each cell's own transform with the parameter name in
  # place of param, not a plausible-looking stand-in: these two cells were
  # exp and log, so that is what it shows.
  expect_match(msg, "MANIFESTMEANS\\[1,1\\] = 'exp\\(mm\\)'")
  expect_match(msg, "MANIFESTMEANS\\[2,1\\] = 'log\\(mm\\)'")
  expect_false(grepl('PARS[row,col]', msg, fixed = TRUE))

  # No transform written anywhere and still two transforms: the DRIFT diagonal
  # default is negative-bounded and the off-diagonal one is identity. This is
  # a plausible thing to write -- "constrain every drift cell to one
  # parameter" -- and it silently produced two different numbers.
  drift <- function(m) suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 2, n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2), DRIFT = m)))
  expect_error(drift(matrix('a', 2, 2)), 'transform it differently')
  expect_error(drift(matrix(c('a', 0, 'a', -0.3), 2, 2)),
    'transform it differently')
  # Across matrices, negative-bounded against positive-bounded.
  expect_error(suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 2, n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    DRIFT = matrix(c('p', 0, 0, -0.3), 2, 2),
    DIFFUSION = matrix(c('p', 0, 0, .4), 2, 2)))), 'transform it differently')
})

test_that("a simple name can be a reference rather than a parameter", {
  # A latent or tdpred name in a cell references that state or predictor, and
  # several matrices may reference the same one -- that is not one parameter
  # appearing twice, so the transform columns of those cells say nothing about
  # each other. The same goes for a PARS label used outside PARS: the cell is
  # a reference, and its own transform is discarded when the reference is
  # rewritten.
  two <- function(...) suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 2, n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), LAMBDA = diag(2), ...)))

  # DRIFT and DIFFUSION default to different transforms, so these would all
  # look like conflicts if a reference were mistaken for a parameter name.
  expect_silent(two(MANIFESTMEANS = matrix(c('eta1', 0), 2, 1),
    CINT = matrix(c('eta1', 0), 2, 1)))
  expect_silent(two(DRIFT = matrix(c('eta1', 0, 0, -0.3), 2, 2),
    DIFFUSION = matrix(c('eta1', 0, 0, .4), 2, 2)))
  expect_silent(two(n.TDpred = 1, TDpredNames = 'TD1',
    DRIFT = matrix(c('TD1', 0, 0, -0.3), 2, 2),
    DIFFUSION = matrix(c('TD1', 0, 0, .4), 2, 2)))
  expect_silent(two(PARS = matrix('p1'),
    DRIFT = matrix(c('p1', 0, 0, -0.3), 2, 2),
    DIFFUSION = matrix(c('p1', 0, 0, .4), 2, 2)))
  expect_silent(two(PARS = matrix('p1|exp(param)'),
    DRIFT = matrix(c('p1', 0, 0, -0.3), 2, 2),
    DIFFUSION = matrix(c('p1', 0, 0, .4), 2, 2)))
})

test_that("declaring the parameter in PARS is accepted, as the message says", {
  # This is the remedy the conflict error recommends, so it has to work: the
  # PARS cell carries the transform and the cells that use the parameter carry
  # only its name.
  m <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('e1', 'e2'),
    LAMBDA = diag(2), PARS = matrix('mm|exp(param)'),
    MANIFESTMEANS = matrix(c('mm', 'mm'), 2, 1)))
  # The using cells become references to the PARS cell, so its transform is
  # the only one that survives.
  rw <- suppressMessages(ctModelStatesAndPARS(m$pars,
    statenames = m$latentNames, tdprednames = m$TDpredNames))
  expect_equal(unique(rw$param[rw$matrix %in% 'MANIFESTMEANS']), 'PARS[1,1]')
  expect_equal(m$pars$transform[m$pars$matrix %in% 'PARS'], 'exp(param)')
})

test_that("two PARS cells sharing a label must still agree on the transform", {
  # Repeating a label across PARS cells is the equality constraint inside
  # PARS, so those cells are co-holders and a disagreement between them is
  # real. An earlier version exempted every PARS label from the check and hid
  # exactly this.
  msg <- tryCatch(suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 2, n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), LAMBDA = diag(2),
    PARS = matrix(c('q|exp(param)', 'q|log(param)'), 2, 1),
    DRIFT = matrix(c('q', 0, 0, -0.3), 2, 2)))),
    error = function(e) conditionMessage(e))
  expect_match(msg, 'transform it differently')
  expect_match(msg, 'PARS\\[1,1\\] as exp\\(param\\)')
  expect_match(msg, 'PARS\\[2,1\\] as log\\(param\\)')
  # And the remedy fits the case: one PARS cell, not two.
  expect_match(msg, 'in one PARS cell rather than two')
})

test_that("a whole constrained matrix reports two cells, not all of them", {
  # DRIFT = matrix('a', 2, 2) is four cells with two transforms. Listing the
  # first two *cells* would show two that agree, which reads as a
  # contradiction of the sentence above the list, so the list is one cell per
  # distinct transform.
  msg <- tryCatch(suppressWarnings(suppressMessages(ctModel(type = 'ct',
    n.latent = 2, n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('eta1', 'eta2'), LAMBDA = diag(2),
    DRIFT = matrix('a', 2, 2)))), error = function(e) conditionMessage(e))
  expect_match(msg, 'DRIFT\\[1,1\\] as ')
  expect_match(msg, 'DRIFT\\[1,2\\] as param')
  # The two diagonal cells agree, so only one of them is listed.
  expect_false(grepl('DRIFT[2,2] as', msg, fixed = TRUE))
  expect_match(msg, '4 cells name it, with 2 different transforms')
})

test_that("sharing a name is still how an equality constraint is written", {
  # The rule is about the resolved transform, not about what was written, so
  # cells whose transforms agree keep working -- which is the common case.
  expect_silent(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    DRIFT = matrix(c('d', 0, 0, 'd'), 2, 2))))
  expect_silent(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    MANIFESTMEANS = matrix(c('mm|exp(param)', 'mm|exp(param)'), 2, 1))))
  # T0MEANS and MANIFESTMEANS share a default, so a name may span them.
  expect_silent(suppressMessages(ctModel(type = 'ct', n.latent = 2,
    n.manifest = 2, manifestNames = c('Y1', 'Y2'),
    latentNames = c('e1', 'e2'), LAMBDA = diag(2),
    T0MEANS = matrix(c('q', 0), 2, 1),
    MANIFESTMEANS = matrix(c('q', 0), 2, 1))))
})
