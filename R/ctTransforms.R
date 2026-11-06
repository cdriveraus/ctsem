# A lower factor of an interval's process noise covariance. A latent no noise
# reaches -- a trend, a fixed initial level -- has a row and column of exact
# zeros there, which `chol()` refuses, so the factor is the rest's, zero there.
.ctGenerateRoot <- function(Q){
  keep <- which(diag(Q) > 0)
  L <- matrix(0, nrow(Q), ncol(Q))
  if(length(keep)) L[keep, keep] <- t(chol(Q[keep, keep, drop = FALSE]))
  L
}
