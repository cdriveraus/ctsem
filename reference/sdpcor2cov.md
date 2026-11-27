# sdcor2cov

Converts a lower triangular matrix with standard deviations on the
diagonal and partial correlations on lower triangle, to a covariance (or
cholesky decomposed covariance)

## Usage

``` r
sdpcor2cov(mat, coronly = FALSE, cholesky = FALSE)
```

## Arguments

- mat:

  input square matrix with std dev on diagonal and lower tri of partial
  correlations.

- coronly:

  if TRUE, ignores everything except the lower triangle and outputs
  correlation.

- cholesky:

  Logical. To return the cholesky decomposition instead of full
  covariance, set to TRUE.

## Value

A covariance matrix, or its lower Cholesky factor with
`cholesky = TRUE`, or a correlation matrix with `coronly = TRUE`.

## Examples

``` r
testmat <- diag(exp(rnorm(5,-3,2)),5) #generate arbitrary std deviations
testmat[row(testmat) > col(testmat)] <- runif((5^2-5)/2, -1, 1) 
print(testmat)
covmat <- sdpcor2cov(testmat) #convert to covariance
cov2cor(covmat) #convert covariance to correlation
```
