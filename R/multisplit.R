#' Variable Selection on Random Sample Splits.
#'
#' Performs repeated variable selection via the lasso on random sample splits.
#'
#' @param x The SNP data matrix, of size \code{nobs x nvar}. Each row represents
#' a subject, each column a SNP.
#' @param y The response vector. It can be continuous or discrete.
#' @param covar NULL or the matrix of covariates one wishes to control for, of
#' size \code{nobs x ncovar}.
#' @param B The number of random splits. Default value is 50.
#'
#' @details The samples are divided into two random splits of approximately
#' equal size. The first subsample is used for variable selection, which is
#' implemented using \pkg{\link{glmnet}}. The first \code{[nobs/6]} variables
#' which enter the lasso path are selected. The procedure is repeated \code{B}
#' times.
#'
#' If one or more covariates are specified, these will be added unpenalized to
#' the regression.
#'
#' @return A data frame with 2 components. A matrix of size \code{B x [nobs/2]}
#' containing the second subsample of each split, and a matrix  of size
#' \code{B x [nobs/6]} containing the selected variables in each split.
#'
#' @examples
#' library(MASS)
#' x <- mvrnorm(60,mu = rep(0,200), Sigma = diag(200))
#' beta <- rep(1,200)
#' beta[c(5,9,3)] <- 3
#' y <- x %*% beta + rnorm(60)
#' res.multisplit <- multisplit(x, y)
#'
#' @references Meinshausen, N., Meier, L. and Buhlmann, P. (2009), P-values for
#' high-dimensional regression, Journal of the American Statistical Association
#' 104, 1671-1681.
#'
#' @name multisplit

multisplit <- function(x, y, covar = NULL, B = 50) {

  # check if the response is binary or continuous
  if ((min(y) == 0) & (max(y) == 1) & (length(unique(y)) == 2)) {
    family = "binomial"
  } else {
    family = "gaussian"
  }

  if (is.numeric(y) == FALSE) {
    y <- as.numeric(y)
  }

  if ((is.matrix(x) == FALSE) || (is.numeric(x) == FALSE)) {
    x <- data.matrix(x)
  }

  if (nrow(x) != length(y)) {
    stop("'x' and 'y' should have an equal number of samples")
  }

  if (nrow(x) != length(y)) {
    stop("'x' and 'y' should have an equal number of samples")
  }

  no.snps <- ncol(x)

  # check if we have covariates
  if (!is.null(covar)) {
    covar <- data.matrix(covar)
    if (length(y) != nrow(covar)) {
      stop("'y' and 'covar' should have an equal number of samples")
    }
    xc <- cbind(covar,x)
    # don't penalize the covariates
    ncovar <- ncol(covar)
    pen <- c(rep(0,ncovar),rep(1,(no.snps - ncovar)))
  } else {
    ncovar <- 0
    pen <- c(rep(1,no.snps))
  }

  # get the lasso coefficients and the remaining samples
  ret <- lapply(1:B, lasso.select, x = x, y = y, ncovar = ncovar, pen = pen,
                  family = family)

  # return the values
  # Initialize the matrix that stores the selected lasso coefficients
  sel.coeff.matrix <- matrix(0, nrow = B, ncol = length(ret[[1]]$sel.coeff))

  # Initialize the matrix that stores the indices of the out.sample
  out.sample.matrix <- matrix(0, nrow = B, ncol = length(ret[[1]]$out.sample))

  for (b in 1:B)
  {
    sel.coeff.matrix[b,] <- ret[[b]]$sel.coeff
    out.sample.matrix[b,] <- ret[[b]]$out.sample
  }

  out.data <- list("out.sample" = out.sample.matrix,
                   "sel.coeff" = sel.coeff.matrix)
  return(out.data)
}
