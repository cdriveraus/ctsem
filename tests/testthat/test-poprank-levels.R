if (!identical(Sys.getenv("NOT_CRAN"), "true")) skip_on_cran()

# A per-level `poprank` under intoverpop='laplace'.
#
# The restriction is applied where the covariance is built -- a loading matrix
# per level in the engine -- rather than by rewriting the model into a basis
# and regressions, which is what the augmented route does and which reaches the
# innermost level only. So these tests assert on the parameter allocation each
# level receives, because that is the thing that differs.

levelframe <- function(nstudy = 6L, nperson = 10L, nburst = 4L) {
  d <- expand.grid(time = c(0, .1, .25, .4), burst_local = seq_len(nburst),
    person = seq_len(nperson), study = seq_len(nstudy))
  d <- d[order(d$study, d$person, d$burst_local, d$time), ]
  d$subject <- (d$study - 1L) * nperson + d$person
  d$burst <- (d$subject - 1L) * nburst + d$burst_local
  set.seed(1)
  d$PA <- stats::rnorm(nrow(d))
  d$NEG <- stats::rnorm(nrow(d))
  d$RNT <- stats::rnorm(nrow(d))
  d
}

levelmodel <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = c("PA", "NEG", "RNT"),
    latentNames = c("PA", "NEG", "RNT"), LAMBDA = diag(3),
    T0MEANS = matrix(paste0("t0_", 1:3), 3, 1), T0VAR = diag(3),
    MANIFESTMEANS = matrix(paste0("mm_", 1:3), 3, 1), MANIFESTVAR = "diag",
    DRIFT = "auto", CINT = matrix(0, 3, 1), DIFFUSION = "auto",
    id = c("burst", "subject", "study"), time = "time", tipredDefault = FALSE)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying_subject <- FALSE
  m$pars$indvarying_study <- FALSE
  diagonal <- m$pars$row == m$pars$col
  m$pars$indvarying[m$pars$matrix == "T0MEANS"] <- TRUE
  wide <- (m$pars$matrix == "MANIFESTMEANS") |
    (m$pars$matrix == "DRIFT" & diagonal) |
    (m$pars$matrix == "DIFFUSION" & diagonal)
  m$pars$indvarying_subject[wide & is.na(m$pars$value)] <- TRUE
  m$pars$indvarying_study[wide & is.na(m$pars$value)] <- TRUE
  m
}

levelspec <- function(d, ...) {
  suppressMessages(ctFit(d, levelmodel(), backend = "julia",
    intoverpop = "laplace", cores = 2, fit = FALSE, ...))
}

bylevel <- function(spec) {
  out <- lapply(spec$laplace$levels, function(lv) data.frame(
    name = lv$name, k = length(lv$re_index),
    rank = if (is.null(lv$rank)) length(lv$re_index) else lv$rank,
    nsd = length(lv$sd_index), ncor = length(lv$cor_index),
    nload = length(lv$load_index), stringsAsFactors = FALSE))
  do.call(rbind, out)
}

test_that("without poprank every level keeps its full covariance", {
  d <- levelframe()
  tab <- bylevel(levelspec(d))
  expect_equal(tab$name, c("burst", "subject", "study"))
  expect_equal(tab$k, c(3L, 9L, 9L))
  expect_equal(tab$rank, tab$k)
  expect_equal(tab$nsd, tab$k)
  expect_equal(tab$ncor, as.integer(tab$k * (tab$k - 1) / 2))
  expect_equal(tab$nload, rep(0L, 3))
})

test_that("a named poprank reduces only the level it names", {
  d <- levelframe()
  full <- levelspec(d)
  reduced <- levelspec(d, poprank = c(study = 2))
  tab <- bylevel(reduced)
  # burst and subject untouched, study a rank-2 loading: 9*2 - 1 = 17.
  expect_equal(tab$rank, c(3L, 9L, 2L))
  expect_equal(tab$nload, c(0L, 0L, 17L))
  expect_equal(tab$nsd, c(3L, 9L, 0L))
  expect_equal(tab$ncor, c(3L, 36L, 0L))
  # Exactly the study level's saving: 9 + 36 scales and correlations replaced
  # by 17 loadings.
  expect_equal(full$laplace$npar - reduced$laplace$npar, 9L + 36L - 17L)
})

test_that("a length-one named poprank is not read as one rank for every level", {
  # `c(study = 2)` is length one and names a level. Reading it as "2
  # everywhere" reduces every level while looking like it did as it was told.
  tab <- bylevel(levelspec(levelframe(), poprank = c(study = 2)))
  expect_equal(tab$rank[tab$name == "burst"], 3L)
  expect_equal(tab$rank[tab$name == "subject"], 9L)
})

test_that("an unnamed single poprank applies to every level", {
  tab <- bylevel(levelspec(levelframe(), poprank = 2))
  expect_equal(tab$rank, c(2L, 2L, 2L))
  # 3*2 - 1 = 5 at the burst level, 9*2 - 1 = 17 at the other two.
  expect_equal(tab$nload, c(5L, 17L, 17L))
  expect_equal(tab$nsd, c(0L, 0L, 0L))
})

test_that("several levels can be reduced at once", {
  tab <- bylevel(levelspec(levelframe(), poprank = c(subject = 3, study = 2)))
  expect_equal(tab$rank, c(3L, 3L, 2L))
  expect_equal(tab$nload, c(0L, 24L, 17L))
})

test_that("a rank at or above the level's own k is the full covariance", {
  tab <- bylevel(levelspec(levelframe(), poprank = c(study = 50)))
  expect_equal(tab$rank[tab$name == "study"], 9L)
  expect_equal(tab$nload[tab$name == "study"], 0L)
})

test_that("poprank refuses what it cannot mean", {
  d <- levelframe()
  expect_error(levelspec(d, poprank = c(nosuchlevel = 2)), "no level called")
  expect_error(levelspec(d, poprank = c(study = 0)), "below 1")
  expect_error(levelspec(d, poprank = c(2, 3)), "must be named")
})

test_that("auto resolves to each level's own mean-affecting count", {
  # MANIFESTMEANS and the DRIFT diagonal reach the observation mean; the
  # DIFFUSION diagonal does not. So the burst level (T0MEANS only) is already
  # full rank and is left alone, while subject and study drop from 9 to 6.
  tab <- bylevel(levelspec(levelframe(), poprank = "auto"))
  expect_equal(tab$rank, c(3L, 6L, 6L))
  expect_equal(tab$nload, c(0L, 9L * 6L - 15L, 9L * 6L - 15L))
  expect_equal(tab$nsd, c(3L, 0L, 0L))
})

test_that("a reduced level's parameters are named as loadings, not scales", {
  d <- levelframe()
  spec <- levelspec(d, poprank = c(study = 2))
  npar <- spec$laplace$npar
  # `.ctBackendRawParameterNames()` reads `fit$model_spec`, and a `fit = FALSE`
  # call returns exactly that object, so this is the same input a fitted object
  # would hand it.
  names <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)
  loadings <- grep("^poploading_", names, value = TRUE)
  expect_equal(length(loadings), 17L)
  expect_true(all(grepl("\\.study$", loadings)))
  # Two dimensions, and nothing at the study level still called a scale or a
  # correlation: reporting a loading as a standard deviation would be a
  # plausible number with the wrong meaning.
  expect_true(any(grepl("_dim1\\.study$", loadings)))
  expect_true(any(grepl("_dim2\\.study$", loadings)))
  expect_false(any(grepl("^popsd_.*\\.study$", names)))
  expect_false(any(grepl("^rawcor_.*\\.study$", names)))
})

test_that("a reduced fit runs and reports a likelihood", {
  d <- levelframe()
  expected <- levelspec(d, poprank = c(subject = 3, study = 2))$laplace$npar
  fit <- suppressMessages(ctFit(d, levelmodel(), backend = "julia",
    intoverpop = "laplace", cores = 2, poprank = c(subject = 3, study = 2),
    optimcontrol = list(estonly = TRUE, maxiter = 5)))
  expect_true(is.finite(as.numeric(fit$estimate$loglik)))
  # The fitted object keeps the spec under `model_spec`; `fit$laplace` is the
  # run's own inner-loop report and carries no `npar`.
  expect_equal(length(fit$estimate$raw), expected)
  expect_equal(fit$model_spec$laplace$levels[[3]]$rank, 2L)
})

test_that("priors cover a reduced level's loadings", {
  # The prior layout is built from each level's scale and correlation blocks,
  # and a reduced level has neither -- so before this was handled, `priors =
  # TRUE` refused the model outright with "the Laplace layout accounts for 45
  # of 156 free parameters". Nothing in a fit without priors would have shown
  # it, which is why this is its own test.
  d <- levelframe()
  fit <- suppressMessages(ctFit(d, levelmodel(), backend = "julia",
    intoverpop = "laplace", cores = 2, priors = TRUE,
    poprank = c(subject = 3, study = 2),
    optimcontrol = list(estonly = TRUE, maxiter = 3)))
  expect_true(is.finite(as.numeric(fit$estimate$loglik)))
})

test_that("a loading prior is scaled so the implied variance does not track rank", {
  d <- levelframe()
  spec <- levelspec(d, poprank = c(study = 2))
  prior <- ctsem:::.ctBackendLaplacePriorSpec(spec$standata, spec$laplace,
    spec$laplace$npar)
  expect_equal(length(prior$index), spec$laplace$npar)
  study <- spec$laplace$levels[[3]]
  expect_equal(unique(prior$scale[match(study$load_index, prior$index)]),
    1 / sqrt(2))
})
