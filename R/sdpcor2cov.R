#' sdcor2cov
#' 
#' Converts a lower triangular matrix with standard deviations on the diagonal and partial correlations on
#' lower triangle, to a covariance (or cholesky decomposed covariance)
#' @param mat input square matrix with std dev on diagonal and lower tri of partial correlations.
#' @param coronly if TRUE, ignores everything except the lower triangle and outputs correlation.
#' @param cholesky Logical. To return the cholesky decomposition instead of full covariance, set to TRUE.
#' @examples
#' testmat <- diag(exp(rnorm(5,-3,2)),5) #generate arbitrary std deviations
#' testmat[row(testmat) > col(testmat)] <- runif((5^2-5)/2, -1, 1) 
#' print(testmat)
#' covmat <- sdpcor2cov(testmat) #convert to covariance
#' cov2cor(covmat) #convert covariance to correlation
#' @export
sdpcor2cov <- function(mat, coronly=FALSE, cholesky=FALSE){ 
  
  ndim = ncol(mat);
  mcholcor=diag(0,ndim);
  mcholcor[1,1]=1;
  
  if(ndim > 1){
    for(coli in 1:ndim){
      for(rowi in coli:ndim){
        if(coli==1 && rowi > 1) mcholcor[rowi,coli] =  mat[rowi,coli]; 
        if(coli > 1){
          if(rowi == coli) mcholcor[rowi,coli] = prod(sqrt(1-mat[rowi,1:(coli-1)]^2));
          if(rowi > coli) mcholcor[rowi,coli] = mat[rowi,coli] * prod(sqrt(1-mat[rowi,1:(coli-1)]^2));
        }
      }
    }
  }
  
  if(!coronly){
    mscale=diag(diag(mat))
    out= mscale %*% mcholcor
  } else out = mcholcor
  if(!cholesky) out = out %*% t(out)
  return(out);
}


constraincorsqrt1 <- function(mat) {
  d <- nrow(mat)
  if (ncol(mat) != d) stop("`mat` must be square.")
  
  # output matrix
  o <- matrix(0, d, d)
  
  # accumulators
  ss <- numeric(d)
  s  <- numeric(d)
  eps <- 1e-5
  
  # 1) compute ss[i] = sum of squares of off‐diag entries in row i (and col i)
  #    and s[i] = sum of off‐diag entries
  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      if (j > i) {
        ss[i] <- ss[i] + mat[j, i]^2
        s[i]  <- s[i] + mat[j, i]
      } else if (j < i) {
        ss[i] <- ss[i] + mat[i, j]^2
        s[i]  <- s[i] + mat[i, j]
      }
    }
    ss[i] <- ss[i] + eps
    s[i]  <- s[i]  + eps
  }
  
  # 2) build each row i of the constrained matrix
  for (i in seq_len(d)) {
    # intermediate scalars
    r1 <- sqrt(ss[i])
    r3 <- abs(s[i]) / r1 - 1
    # log1p_exp(x) ≈ log(1 + exp(x))
    r4 <- sqrt(log1p(exp(2 * (abs(s[i]) - s[i] - 1) - 4)))
    tmp <- (r4 * r3 + 1) * r4 + 1
    r <- sqrt(ss[i] + tmp)
    
    # fill off‐diagonals in row i
    for (j in seq_len(d)) {
      if (j > i) {
        o[i, j] <- mat[j, i] / r
      } else if (j < i) {
        o[i, j] <- mat[i, j] / r
      }
    }
    
    # set diagonal so that row norm ≈ 1
    # guard against tiny negatives under the sqrt
    diag_val <- 1 - sum(o[i, ]^2) + eps
    o[i, i] <- sqrt(max(0, diag_val))
  }
  
  o
}

constraincorsqrt2 <- function(mat) { #symmetric input form
  d <- nrow(mat)
  if (ncol(mat) != d) stop("`mat` must be square.")
  
  # output matrix
  o <- matrix(0, d, d)
  
  # accumulators
  ss <- numeric(d)
  s  <- numeric(d)
  eps <- 1e-5
  
  # 1) compute ss[i] = sum of squares of off‐diag entries in row i (and col i)
  #    and s[i] = sum of off‐diag entries
  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      # if (j > i) {
        ss[i] <- ss[i] + mat[j, i]^2
        s[i]  <- s[i] + mat[j, i]
      # } else if (j < i) {
      #   ss[i] <- ss[i] + mat[i, j]^2
      #   s[i]  <- s[i] + mat[i, j]
      # }
    }
    ss[i] <- ss[i] + eps
    s[i]  <- s[i]  + eps
  }
  
  # 2) build each row i of the constrained matrix
  for (i in seq_len(d)) {
    # intermediate scalars
    r1 <- sqrt(ss[i])
    r3 <- abs(s[i]) / r1 - 1
    # log1p_exp(x) ≈ log(1 + exp(x))
    r4 <- sqrt(log1p(exp(2 * (abs(s[i]) - s[i] - 1) - 4)))
    tmp <- (r4 * r3 + 1) * r4 + 1
    r <- sqrt(ss[i] + tmp)
    
    # fill off‐diagonals in row i
    for (j in seq_len(d)) {
      # if (j > i) {
      #   o[i, j] <- mat[j, i] / r
      # } else if (j < i) {
        o[i, j] <- mat[i, j] / r
      # }
    }
    
    # set diagonal so that row norm ≈ 1
    # guard against tiny negatives under the sqrt
    diag_val <- 1 - sum(o[i, ]^2) + eps
    o[i, i] <- sqrt(max(0, diag_val))
  }
  
  o
}


