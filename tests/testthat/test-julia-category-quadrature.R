test_that('ordinal category probabilities match independent integration', {
  skip_on_cran()
  skip_without_julia()
  ctJuliaSetup()

  # The engine's categorical measurement update is an adaptive Gauss-Hermite
  # rule whose placement comes from a Newton solve for the posterior mode. When
  # that solve misses, the rule is still a proper quadrature -- it is merely
  # centred where the posterior has no mass, so the category's probability comes
  # back wrong and nothing errors. This asserts the value against an
  # independently derived answer rather than against another run of the same
  # code, which is the only thing that catches it: a parity test between two
  # callers of this kernel agrees with itself while both are wrong.
  #
  # The regression it guards: the Newton step is `gradient / curvature`, and an
  # ordinal category the predicted state makes improbable has a log likelihood
  # asymptotically linear in the state -- score near +-1, information decaying
  # to zero. The curvature left is then the prior's 1/s^2 and the step is
  # +-s^2, which oscillates instead of converging. Measured at s = 10 before
  # `_mode_start` and the step cap: 27 log units out on an interior category,
  # against 3e-06 for the categories the solve did reach. More quadrature nodes
  # barely move it, which is what separates this from an under-resolved
  # integral.
  JuliaConnectoR::juliaEval('
    function _ctsem_test_logp(etabar, s, y, tau, n)
      M = ContinuousTimeSEM
      nd, wt = M._gauss_hermite(Int(n))
      lz, _, _ = M._binary_moments(etabar, s, Int(y), nd, wt, tau,
        M.CTSEM_OBS_ORDINAL)
      lz
    end')

  # P(y = k) = E_eta[ sigma(tau_k - eta) - sigma(tau_{k-1} - eta) ], by adaptive
  # quadrature in R: a different algorithm, in a different language, from the
  # cumulative logit as written rather than as the engine arranges it.
  reference <- function(etabar, s, y, tau) {
    k <- length(tau) + 1L
    upper <- if (y == k) Inf else tau[y]
    lower <- if (y == 1L) -Inf else tau[y - 1L]
    integrand <- function(e) {
      ((if (is.infinite(upper)) 1 else stats::plogis(upper - e)) -
          (if (is.infinite(lower)) 0 else stats::plogis(lower - e))) *
        stats::dnorm(e, etabar, s)
    }
    stats::integrate(integrand, etabar - 40 * s, etabar + 40 * s,
      rel.tol = 1e-12, subdivisions = 4000L)$value
  }

  tau <- seq(-3, 3, length.out = 8)          # nine categories
  ncat <- length(tau) + 1L

  # s is the predicted sd of the linear predictor. Small s was always fine; the
  # failure needs a diffuse predicted state, which an optimiser reaches whenever
  # it tries a large DIFFUSION or T0VAR, so these are ordinary values rather
  # than a corner.
  for (s in c(2, 5, 10)) {
    for (etabar in c(-3, 0, 3)) {
      got <- vapply(seq_len(ncat), function(y)
        JuliaConnectoR::juliaCall('_ctsem_test_logp', etabar, s, y, tau, 21),
        numeric(1))
      want <- vapply(seq_len(ncat), function(y)
        log(reference(etabar, s, y, tau)), numeric(1))
      # Elementwise, so a failure names the category that moved. 1e-3 is far
      # above the rule's own error here (worst measured 3.8e-04 at s = 10) and
      # far below the defect (27).
      expect_equal(got, want, tolerance = 1e-3,
        info = paste0('s = ', s, ', etabar = ', etabar))
    }
  }
})

test_that('the ordinal mode solve lands in the observed category band', {
  skip_on_cran()
  skip_without_julia()
  ctJuliaSetup()

  # The mechanism, asserted directly rather than through its consequence. With
  # a diffuse prior the posterior mode of an interior category sits inside that
  # category's own band; the divergent solve returned +-s^2 instead, and the
  # sign of it flipped with the iteration count.
  JuliaConnectoR::juliaEval('
    function _ctsem_test_mode(etabar, s, y, tau)
      M = ContinuousTimeSEM
      mo, _ = M._binary_mode(etabar, s * s, Int(y), tau, M.CTSEM_OBS_ORDINAL)
      mo
    end')

  tau <- seq(-3, 3, length.out = 8)
  for (s in c(5, 10, 25)) {
    for (y in 2:8) {                      # interior categories only
      offset <- JuliaConnectoR::juliaCall('_ctsem_test_mode', 0, s, y, tau)
      # A margin either side of the band: with a finite prior the mode may sit
      # slightly outside a narrow band, but not by the units a runaway step
      # produces.
      expect_gt(offset, tau[y - 1L] - 2)
      expect_lt(offset, tau[y] + 2)
    }
  }
})

test_that('every categorical kind integrates to an independent reference', {
  skip_on_cran()
  skip_without_julia()
  ctJuliaSetup()

  # The ordinal test above covers the kind the defect was found in. This one
  # covers the other three, because the mode solve they share is the same code
  # and the failure it had was not ordinal-specific: a bounded score with a
  # vanishing information produces a step of the prior variance for any
  # likelihood whose log flattens, which is all of them somewhere.
  #
  # Tolerances are per kind and per predictor sd because the remaining error is
  # the quadrature's own and it differs between them. They are set from
  # measurement, a little above what the engine currently achieves, so that a
  # regression shows up as a failure rather than as a slow drift. Where they
  # look loose -- censored at a predictor sd of 10 -- that is a real limit of a
  # 21 node rule against a near-kinked integrand, recorded here rather than
  # hidden.
  JuliaConnectoR::juliaEval('
    function _ctsem_test_kind(etabar, s, y, tau, kind)
      M = ContinuousTimeSEM
      nd, wt = M._gauss_hermite(M._CTSEM_BINARY_NODES[])
      lz, _, _ = M._binary_moments(etabar, s, y, nd, wt, tau, Int(kind))
      lz
    end')
  logp <- function(etabar, s, y, tau, kind)
    JuliaConnectoR::juliaCall('_ctsem_test_kind', etabar, s, y, tau, kind)

  ## ---- binary: E[sigma(eta)] by adaptive integration
  for (case in list(list(s = 2, tol = 1e-5), list(s = 5, tol = 5e-3),
    list(s = 10, tol = 2e-2))) {
    for (etabar in c(-4, 0, 4)) for (y in c(0, 1)) {
      want <- log(stats::integrate(function(x)
        (if (y > 0.5) stats::plogis(x) else stats::plogis(-x)) *
          stats::dnorm(x, etabar, case$s),
        etabar - 40 * case$s, etabar + 40 * case$s,
        rel.tol = 1e-12, subdivisions = 4000L)$value)
      expect_equal(logp(etabar, case$s, y, numeric(0), 1), want,
        tolerance = case$tol,
        info = paste0('binary s = ', case$s, ', etabar = ', etabar, ', y = ', y))
    }
  }

  ## ---- censored: Gaussian on a Gaussian, so the engine solves it in closed
  ## form and the reference is closed form too. The tolerance is at machine
  ## precision on purpose: it asserts that the analytic path is the one being
  ## taken, and would fail immediately if a censored row ever fell back to the
  ## quadrature, which was 7.4e-02 out at the last of these.
  limits <- c(-2, 2, 1)                      # lower, upper, measurement sd
  for (s in c(0.5, 2, 5, 10, 20)) {
    total <- sqrt(1 + s^2)
    for (etabar in c(-6, 0, 6)) for (y in c(-2, 0, 2)) {
      want <- if (y == -2) stats::pnorm(-2, etabar, total, log.p = TRUE) else
        if (y == 2) stats::pnorm(2, etabar, total, lower.tail = FALSE,
          log.p = TRUE) else stats::dnorm(y, etabar, total, log = TRUE)
      expect_equal(logp(etabar, s, y, limits, 4), want, tolerance = 1e-12,
        info = paste0('censored s = ', s, ', etabar = ', etabar, ', y = ', y))
    }
  }

  ## ---- count: Poisson with a lognormal rate. There is no closed form, so
  ## the reference is computed two ways -- adaptively and on a fine fixed grid
  ## -- and a cell where they disagree is skipped rather than compared. A
  ## deeply improbable count puts the integrand many prior sd out, where
  ## adaptive integration can miss the contributing region entirely and return
  ## a confident wrong answer; comparing the engine against that would be
  ## testing R's quadrature, not ctsem's.
  countref <- function(etabar, s, y) {
    f <- function(x) stats::dpois(y, exp(x)) * stats::dnorm(x, etabar, s)
    a <- try(stats::integrate(f, etabar - 40 * s, etabar + 40 * s,
      rel.tol = 1e-13, subdivisions = 6000L)$value, silent = TRUE)
    if (inherits(a, 'try-error') || !is.finite(a) || a <= 0) return(NA_real_)
    grid <- seq(etabar - 25 * s, etabar + 25 * s, length.out = 100001)
    b <- sum(f(grid)) * (grid[2] - grid[1])
    if (!is.finite(b) || b <= 0) return(NA_real_)
    if (abs(log(a) - log(b)) > 1e-7) return(NA_real_)
    log(a)
  }
  checked <- 0L
  for (case in list(list(s = 0.69, tol = 1e-6), list(s = 1.2, tol = 1e-5),
    list(s = 2, tol = 1e-3))) {
    for (etabar in c(-3, 0, 1.11, 3)) for (y in c(0, 1, 3, 20)) {
      want <- countref(etabar, case$s, y)
      if (is.na(want)) next
      checked <- checked + 1L
      expect_equal(logp(etabar, case$s, y, numeric(0), 3), want,
        tolerance = case$tol,
        info = paste0('count s = ', case$s, ', etabar = ', etabar, ', y = ', y))
    }
  }
  # The skipping above is a real hazard: were the reference to become unusable
  # everywhere, this block would pass by testing nothing at all.
  expect_gt(checked, 25L)
})

test_that('the mode solve reaches a stationary point for every kind', {
  skip_on_cran()
  skip_without_julia()
  ctJuliaSetup()

  # The property the line search is there to deliver, asserted directly: at the
  # returned offset the gradient of `log N(eta; etabar, s^2) + log P(y | eta)`
  # is zero to the tolerance the solve stops at. This is what a fixed iteration
  # count could not promise and what neither a per-kind step cap nor a longer
  # budget delivered -- the iteration it replaced oscillated, so it could sit a
  # long way from stationary with no sign of it.
  JuliaConnectoR::juliaEval('
    function _ctsem_test_stationary(etabar, s, y, tau, kind)
      M = ContinuousTimeSEM
      k = Int(kind)
      mo, curv = M._binary_mode(etabar, s * s, y, tau, k)
      score, _ = M._category_score(etabar + mo, y, tau, k)
      gradient = -mo / (s * s) + score
      # The objective still available, which is what the solve stops on and the
      # only form of this that is comparable across kinds and across s: a raw
      # gradient bound asks for a number of digits that depends on the
      # curvature, which is the mistake `_laplace_inner_tolerance` is written
      # to avoid.
      abs(gradient * gradient / (2 * curv))
    end')
  remaining <- function(etabar, s, y, tau, kind)
    JuliaConnectoR::juliaCall('_ctsem_test_stationary', etabar, s, y, tau, kind)

  tau <- seq(-3, 3, length.out = 8)
  # Three orders above the solve's own tolerance, so this asserts that it
  # stopped because it had converged rather than because it ran out of
  # iterations, without pinning the exact stopping rule.
  for (s in c(0.5, 2, 5, 10, 25)) {
    for (etabar in c(-6, -2, 0, 2, 6)) {
      for (y in 1:9) expect_lt(remaining(etabar, s, y, tau, 2), 1e-9)
      for (y in c(0, 1)) expect_lt(remaining(etabar, s, y, numeric(0), 1), 1e-9)
      for (y in c(0, 1, 20))
        expect_lt(remaining(etabar, s, y, numeric(0), 3), 1e-9)
      for (y in c(-2, 0, 2))
        expect_lt(remaining(etabar, s, y, c(-2, 2, 1), 4), 1e-9)
    }
  }
})

test_that('binary asymptotes reduce to plain binary and match direct integration', {
  skip_on_cran()
  skip_without_julia()
  ctJuliaSetup()

  # A three or four parameter logistic likelihood is not log-concave for a
  # correct response, so its scalar posterior can have two modes -- measured at
  # a prior mean below about -4.3 with c = 0.2. The quadrature path assumes one
  # mode and would centre on whichever it found. `_asymptote_moments` avoids
  # the question rather than answering it: `P(y|eta) = A + B q(y|eta)` is a
  # mixture of the prior and the plain logistic posterior, both unimodal, so
  # the components are integrated apart and their moments combined.
  #
  # Two things have to hold. Without asymptotes nothing may change, and with
  # them the mixture must agree with integrating the lump -- including where
  # the lump is bimodal, which is the case the whole construction is for.
  JuliaConnectoR::juliaEval('
    function _ctsem_test_bin(etabar, s, y, extras)
      M = ContinuousTimeSEM
      nd, wt = M._gauss_hermite(M._CTSEM_BINARY_NODES[])
      lz, off, v = M._binary_moments(etabar, s, y, nd, wt, extras,
        M.CTSEM_OBS_BINARY)
      [lz, off, v]
    end')
  moments <- function(etabar, s, y, extras)
    JuliaConnectoR::juliaCall('_ctsem_test_bin', etabar, s, y, extras)

  for (etabar in c(-4, -1, 0, 2, 5)) for (s in c(0.3, 1, 3, 10))
    for (y in c(0, 1)) {
      expect_equal(moments(etabar, s, y, c(0, 1)),
        moments(etabar, s, y, numeric(0)), tolerance = 1e-10,
        info = paste0('reduction at etabar = ', etabar, ', s = ', s))
    }

  # P(y=1) = c + (d-c) plogis(eta), integrated directly on a grid fine enough
  # to resolve it and wide enough to hold the prior.
  direct <- function(etabar, s, y, cc, dd) {
    e <- seq(etabar - 16 * s, etabar + 16 * s, length.out = 200001)
    w <- e[2] - e[1]
    p <- cc + (dd - cc) * stats::plogis(e)
    f <- (if (y == 1) p else 1 - p) * stats::dnorm(e, etabar, s)
    Z <- sum(f) * w
    mu <- sum(e * f) * w / Z
    c(log(Z), mu - etabar, sum(e * e * f) * w / Z - mu^2)
  }
  # The starred rows are prior means where the posterior has two modes.
  cases <- list(
    list(-4.7, 1, 1, 0.20, 1.00),      # bimodal
    list(-12,  2, 1, 0.20, 1.00),      # bimodal, modes ten units apart
    list(0,    1, 1, 0.20, 1.00),
    list(0,    1, 0, 0.20, 1.00),
    list(2,    2, 1, 0.25, 0.90),      # 4PL
    list(-6,   2, 0, 0.15, 0.85),      # 4PL, incorrect response
    list(0,  0.5, 1, 0.35, 1.00))
  for (cs in cases) {
    got <- moments(cs[[1]], cs[[2]], cs[[3]], c(cs[[4]], cs[[5]]))
    want <- direct(cs[[1]], cs[[2]], cs[[3]], cs[[4]], cs[[5]])
    # Loose enough for the plain binary rule's own error at these predictor
    # sds, which the mixture inherits and does not add to.
    expect_equal(got, want, tolerance = 1e-4,
      info = paste0('etabar = ', cs[[1]], ', s = ', cs[[2]], ', y = ', cs[[3]],
        ', c = ', cs[[4]], ', d = ', cs[[5]]))
  }
})

test_that('the asymptote score is the derivative of the asymptote likelihood', {
  skip_on_cran()
  skip_without_julia()
  ctJuliaSetup()

  # The pair has to describe the same function: the adjoint's
  # degenerate-variance branch takes the score from here while the forward pass
  # takes the value from there, and a mismatch between them is a wrong gradient
  # with a right likelihood.
  #
  # The information is deliberately not clamped for this kind. It is genuinely
  # negative where the likelihood is convex, which is what non-log-concavity
  # means, and reporting that as a floor would be reporting a convex region as
  # a flat one.
  JuliaConnectoR::juliaEval('
    function _ctsem_test_ll(eta, y, extras)
      ContinuousTimeSEM._category_loglikelihood(eta, y, extras, 1) end
    function _ctsem_test_sc(eta, y, extras)
      s, i = ContinuousTimeSEM._category_score(eta, y, extras, 1); [s, i] end')
  ll <- function(e, y, ex) JuliaConnectoR::juliaCall('_ctsem_test_ll', e, y, ex)
  sc <- function(e, y, ex) JuliaConnectoR::juliaCall('_ctsem_test_sc', e, y, ex)

  h <- 1e-4
  sawnegative <- FALSE
  for (cc in c(0, 0.15, 0.3)) for (dd in c(0.85, 1))
    for (e in c(-3, -1, 0, 1, 3)) for (y in c(0, 1)) {
      ex <- c(cc, dd)
      got <- sc(e, y, ex)
      expect_equal(got[1], (ll(e + h, y, ex) - ll(e - h, y, ex)) / (2 * h),
        tolerance = 1e-5, info = paste0('score c=', cc, ' d=', dd, ' eta=', e))
      expect_equal(got[2],
        -(ll(e + h, y, ex) - 2 * ll(e, y, ex) + ll(e - h, y, ex)) / h^2,
        tolerance = 1e-4,
        info = paste0('information c=', cc, ' d=', dd, ' eta=', e))
      if (got[2] < 0) sawnegative <- TRUE
    }
  # If this stops being true the grid has drifted away from the convex region
  # and the test is no longer exercising what it was written for.
  expect_true(sawnegative)

  for (e in c(-6, -2, 0, 2, 6)) for (y in c(0, 1))
    expect_equal(ll(e, y, c(0, 1)), ll(e, y, numeric(0)), tolerance = 1e-14)
})

test_that('the gradient of a binary model with asymptotes matches finite differences', {
  skip_on_cran()
  skip_without_julia()

  # The forward pass and the reverse pass are separate code, and a reverse pass
  # can be wrong while every likelihood in the package is right. That is what
  # happened here: the adjoint mapped an asymptote's cotangent back onto its
  # parameter cell with the rule the ordinal thresholds use -- a reverse
  # cumulative sum, correct when threshold k is a sum of gaps and wrong when
  # `d = c + (1-c)g`. Nothing in the likelihood tests could see it. The fit
  # stopped at a gradient norm of 6536, where a converged one is 1e-3, and the
  # diagnostics then reported the guessing parameter as unidentified with
  # negative curvature -- a description of a point that was not a mode, which
  # reads exactly like the identification problem the three parameter model is
  # known for. With the mapping corrected the same fit converges, the Hessian
  # is negative definite, and the guessing parameter comes back at 0.199
  # (sd 0.008) against a generating 0.20.
  #
  # So: every kind whose parameters reach the engine through the extras slot
  # needs a finite difference check on the *gradient*, not only on the value.
  nit <- 5
  nm <- paste0('y', seq_len(nit))
  loadings <- rep(1.3, nit)
  difficulty <- seq(-1.2, 1.2, length.out = nit)
  gen <- ctModel(type = 'ct', n.latent = 1, n.manifest = nit,
    manifestNames = nm, latentNames = 'eta', manifesttype = rep(1L, nit),
    asymptotes = rep(1L, nit), LAMBDA = matrix(loadings, nit, 1),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(1.3), T0VAR = matrix(1.3),
    T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(-loadings * difficulty, nit, 1),
    MANIFESTVAR = diag(0, nit), Tpoints = 5)
  fixed <- gen$pars$matrix %in% 'THRESHOLDS' & gen$pars$col == 1
  gen$pars$value[fixed] <- 0.2
  gen$pars$param[fixed] <- NA
  gen$pars$transform[fixed] <- NA
  set.seed(5)
  d <- data.frame(ctGenerate(gen, n.subjects = 30, Tpoints = 5, dtmean = 0.8,
    backend = 'julia'))

  JuliaConnectoR::juliaEval('
    function _ctsem_test_value(obj, x)
      ContinuousTimeSEM.ctsem_adjoint_gradient(obj, x).value end
    function _ctsem_test_grad(obj, x)
      collect(ContinuousTimeSEM.ctsem_adjoint_gradient(obj, x).gradient) end')

  for (asym in c(0L, 1L, 2L)) {
    model <- ctModel(type = 'ct', n.latent = 1, n.manifest = nit,
      manifestNames = nm, latentNames = 'eta', manifesttype = rep(1L, nit),
      asymptotes = if (asym == 0L) NULL else rep(asym, nit),
      LAMBDA = matrix(c(1, paste0('a_', nm[-1])), nit, 1),
      CINT = matrix(0), T0MEANS = matrix(0), MANIFESTVAR = diag(0, nit))
    model$pars$indvarying <- FALSE
    # One shared asymptote, which is both the usual way to make it estimable
    # and the case where a wrong cotangent is summed over items rather than
    # cancelling.
    if (asym >= 1L) model$pars$param[model$pars$matrix %in% 'THRESHOLDS' &
        model$pars$col == 1 & !is.na(model$pars$param)] <- 'guess'
    if (asym >= 2L) model$pars$param[model$pars$matrix %in% 'THRESHOLDS' &
        model$pars$col == 2 & !is.na(model$pars$param)] <- 'upper'

    fit <- ctFit(d[, c('id', 'time', nm)], model, backend = 'julia', cores = 1,
      optimcontrol = list(maxiter = 1), priors = TRUE)
    objective <- .ctJuliaObjective(fit)
    # Away from the optimum on purpose: at a mode every gradient is near zero
    # and agrees with anything.
    set.seed(1)
    at <- as.numeric(fit$estimate$raw) +
      stats::rnorm(length(fit$estimate$raw), 0, 0.3)
    analytic <- JuliaConnectoR::juliaCall('_ctsem_test_grad', objective, at)
    step <- 1e-5
    numeric <- vapply(seq_along(at), function(i) {
      up <- at; up[i] <- up[i] + step
      down <- at; down[i] <- down[i] - step
      (JuliaConnectoR::juliaCall('_ctsem_test_value', objective, up) -
          JuliaConnectoR::juliaCall('_ctsem_test_value', objective, down)) /
        (2 * step)
    }, numeric(1))
    expect_equal(analytic, numeric, tolerance = 1e-5,
      info = paste0('asymptotes = ', asym))
  }
})
