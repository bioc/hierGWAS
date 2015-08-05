#' R2 computation
#'
#' Calculates the R2 of a cluster of SNPs.
#'
#' @param x The input matrix, of dimension \code{nobs x nvar}. Each row
#' represents a subject, each column a SNP.
#' @param y The response vector. It can be continuous or discrete.
#' @param res.multisplit The output of \code{multisplit}.
#' @param covar \code{NULL} or the matrix of covariates one wishes to control
#' for, of size \code{nobs x ncovar}.
#' @param SNP_index \code{NULL} or the index vector of the cluster of SNPs whose R2 will be
#' computed. See the 'Details' section.
#'
#' @details The R2 of a cluster of SNPs is computed on the second half-samples.
#' The cluster members, are intersected with the SNPs selected by the lasso, and
#' the R2 of this model is calculated. Thus, we obtain B R2 values. Finally, the
#' mean of these values is taken. If the value of \code{SNP_index} is \code{NULL}, the R2 of the full
#' model with all the SNPs will be computed.
#'
#' @return The R2 value of the SNP cluster
#'
#' @examples
#' library(MASS)
#' x <- mvrnorm(60,mu = rep(0,60), Sigma = diag(60))
#' beta <- rep(0,60)
#' beta[c(5,9,3)] <- 1
#' y <- x %*% beta + rnorm(60)
#' SNP_index <- c(5,9,3)
#' res.multisplit <- multisplit(x, y)
#' r2 <- compute.r2(x, y, res.multisplit, SNP_index = SNP_index)
#'
#' @references Buzdugan, L. et al. (2015), Assessing statistical significance in
#' predictive genome-wide association studies. (unpublished)
#'
#' @name compute.r2
#' @aliases compute.r2

compute.r2 <- function(x, y, res.multisplit, covar = NULL, SNP_index = NULL) {

  if (is.numeric(y) == FALSE) {
    y <- as.numeric(y)
  }
  if ((is.matrix(x) == FALSE) || (is.numeric(x) == FALSE)) {
    x <- data.matrix(x)
  }
  if (nrow(x) != length(y)) {
    stop("'x' and 'y' should have an equal number of samples")
  }
  if (!is.null(covar)) {
    if (length(y) != nrow(covar)) {
      stop("'y' and 'covar' should have an equal number of samples")
    }
  }
  # get the samples and coefficients
  out.sample.matrix <- res.multisplit$out.sample
  lasso.coeff.matrix <- res.multisplit$sel.coeff
  B <- nrow(lasso.coeff.matrix)
  ret.vector <- vector("numeric", length = B)
  # compute the R2 value for each split
  for (b in 1:B) {
    scl <- lasso.coeff.matrix[b,]
    scl <- scl[!is.na(scl)]
    if (is.null(SNP_index)) {
        common.SNP_index <- scl
    } else {
        common.SNP_index <- intersect(SNP_index,scl)   
    }
    if (length(common.SNP_index) == 0) {
      ret.vector[b] <- 0
    } else {
      x.b <- x[out.sample.matrix[b, ], common.SNP_index]
      y.b <- y[ out.sample.matrix[b, ]]
      if (is.null(covar)) {
        ret.vector[b] <- return.r2(x.b, y.b)
      } else {
        covar <- data.matrix(covar)
        ret.vector[b] <- return.r2(x.b, y.b,
                                   covar = covar[out.sample.matrix[b, ],])
      }
      if (ret.vector[b] < 0) {
        ret.vector[b] <- 0
      }
    }
  }
  # compute the final R2
  adjusted.ret <- mean(ret.vector)
  return(adjusted.ret)
}


