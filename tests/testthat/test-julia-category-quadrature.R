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
