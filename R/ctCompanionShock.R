# What "a one unit change in process c" brings with it -------------------------
#
# Every version of the model implied regression is the same calculation --
# dtDRIFT(u) %*% C -- differing only in the companion matrix C, whose column c
# says what else moves when process c moves by one unit. Writing them this way
# makes the choice explicit and keeps one code path:
#
#   'experimental'  C = I. Nothing else moves, because the change was imposed.
#                   This is the partial regression E[x(t+u) | x(t)]: the whole
#                   state is conditioned on, so nothing is left to average over
#                   and the answer is the transition matrix itself. It is the
#                   only variant that is a property of the dynamics alone.
#
#   'observational' C = Sigma diag(Sigma)^-1, Sigma the stationary state
#                   covariance (asymDIFFUSIONcov). We *observe* process c one
#                   unit above its expectation and do not hold the others
#                   fixed, so they come along by E[x_r | x_c = 1] =
#                   Sigma_rc / Sigma_cc. This is the simple (marginal)
#                   regression -- what regressing the observed data on x_c
#                   alone actually recovers.
#
#   'shock'         C = Q diag(Q)^-1, Q the diffusion (innovation) covariance.
#                   A system noise innovation of one unit arrives in c, and the
#                   companion innovations are E[w_r | w_c = 1] = Q_rc / Q_cc.
#                   An impulse response rather than a regression: it asks what
#                   one hypothetical shock does, not what the data would show.
#
#   'orthogonal'    C from the Cholesky factor of Q, column-normalised so the
#                   own effect is one unit. The SVAR convention: shocks are made
#                   uncorrelated by attributing shared variance to whichever
#                   process comes first. That makes it depend on the order of
#                   the latents, which is a modelling assumption, not a
#                   property of the fit -- so it is offered but never a default.
#
# Two things worth noticing about these matrices.
#
# The companion matrices are *asymmetric*, and must be. E[w_2 | w_1 = 1] and
# E[w_1 | w_2 = 1] are different numbers unless the two variances match: a one
# unit move in a wide-spread process implies more movement in a narrow one than
# the reverse. A correlation matrix is symmetric, so it cannot be any of these,
# which is why using one was wrong whenever the processes differed in scale.
#
# 'observational' and 'shock' use *different covariance matrices* because they
# ask different questions. How states covary in the long run is not how one
# innovation relates to another; they coincide only under isotropic decay with
# no cross effects.
#
# Finally, two of these are quantities that already have names. With C = I the
# result is the discrete time autoregression/cross-lagged matrix. With
# C = Sigma diag(Sigma)^-1 and the standardisation applied, the result is the
# model implied cross-correlation function Cor(x_r(t+u), x_c(t)) -- the
# counterpart to what ctACF() computes from the data, and the latent
# correlation matrix at t = 0. Nothing further is missing from the family;
# `cov=TRUE` changes the output type rather than the interpretation.

.ctCompanionTypes <- c("experimental", "observational", "shock", "orthogonal")

# `observational` accepts the historical logical as well as a name.
.ctCompanionType <- function(observational) {
  if (is.logical(observational)) {
    return(if (isTRUE(observational)) "observational" else "experimental")
  }
  match.arg(as.character(observational)[1L], .ctCompanionTypes)
}

# Regression of every process on process c, from a covariance matrix:
# column c is cov[, c] / cov[c, c], so the diagonal is one by construction.
#
# A process with no variance of the relevant kind -- a deterministic latent, a
# process with no diffusion -- gets no companions rather than a division by
# zero, which is the right answer as well as a safe one.
.ctCovarianceRegression <- function(covariance, nlatent) {
  covariance <- as.matrix(covariance)[seq_len(nlatent), seq_len(nlatent), drop = FALSE]
  variance <- diag(covariance)
  variance[!is.finite(variance) | variance <= 0] <- Inf
  out <- covariance %*% diag(1 / variance, nlatent)
  out[!is.finite(out)] <- 0
  diag(out) <- 1
  out
}

#' The companion matrix for one interpretation of a unit change
#'
#' @param type One of .ctCompanionTypes.
#' @param diffusion DIFFUSIONcov, the innovation covariance.
#' @param asymdiffusion asymDIFFUSIONcov, the stationary state covariance.
#' @return nlatent by nlatent matrix, or NULL when the type cannot be formed
#'   from the matrices supplied (a non-stationary linearisation, say).
#' @noRd
.ctCompanionMatrix <- function(type, diffusion, asymdiffusion, nlatent) {
  # Validated here as well as in .ctCompanionType: the branches below end in a
  # fall-through, so an unrecognised name would silently return the last one.
  type <- match.arg(type, .ctCompanionTypes)
  if (identical(type, "experimental")) return(diag(nlatent))

  if (identical(type, "observational")) {
    variance <- diag(as.matrix(asymdiffusion)[seq_len(nlatent), seq_len(nlatent), drop = FALSE])
    # No stationary covariance, no observational interpretation: there is no
    # distribution of states for the companions to be an expectation over.
    if (any(!is.finite(variance)) || any(variance < 0)) return(NULL)
    return(.ctCovarianceRegression(asymdiffusion, nlatent))
  }

  if (identical(type, "shock")) return(.ctCovarianceRegression(diffusion, nlatent))

  # 'orthogonal': Q = L L', so a unit shock to the jth orthogonalised component
  # moves the system by L[, j]. Normalised by its own entry so the own effect is
  # one unit, matching every other column here.
  Q <- as.matrix(diffusion)[seq_len(nlatent), seq_len(nlatent), drop = FALSE]
  L <- try(t(chol(Q + diag(1e-10, nlatent))), silent = TRUE)
  if (inherits(L, "try-error")) return(NULL)
  own <- diag(L)
  own[!is.finite(own) | own == 0] <- Inf
  out <- L %*% diag(1 / own, nlatent)
  out[!is.finite(out)] <- 0
  out
}

# One line naming what was computed, for the help and for messages.
.ctCompanionDescription <- function(type) {
  switch(type,
    experimental = paste0("independent unit impulses (partial regression): ",
      "the other processes are held where they were"),
    observational = paste0("an observed unit change (simple regression): the ",
      "other processes move with it by Sigma_rc/Sigma_cc"),
    shock = paste0("one correlated system noise innovation: the other ",
      "processes get companion shocks of Q_rc/Q_cc"),
    orthogonal = paste0("orthogonalised shocks from the Cholesky factor of ",
      "DIFFUSION -- these depend on the order of the latent processes"))
}
