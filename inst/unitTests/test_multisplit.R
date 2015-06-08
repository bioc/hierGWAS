test_multisplit <- function() {
  library(MASS)
  library(RUnit)
  x <- mvrnorm(50,mu = rep(0,100), Sigma = diag(100))
  beta <- rep(0,100)
  beta[1:5] <- 100
  y <- x%*%beta + rnorm(50)
  covar <- mvrnorm(5,mu = rep(0,50), Sigma = diag(50))

  labels <- seq(1,10)

  out <- multisplit(x,y)

  checkEquals(ncol(out$out.sample),25)
  checkEquals(nrow(out$out.sample),50)
  checkEquals(ncol(out$sel.coeff),floor(50/6))
  checkEquals(nrow(out$sel.coeff),50)
}

