#' Hierarchical Testing of SNPs
#'
#' Performs hierarchical testing of SNPs.
#'
#' @param x The input matrix, of dimension \code{nobs x nvar}. Each row
#' represents a subject, each column a SNP.
#' @param y The response vector. It can be continuous or discrete.
#' @param dendr The cluster tree obtained by hierchically clustering the SNPs
#' using \code{cluster.snp}.
#' @param res.multisplit The output of \code{multisplit}.
#' @param covar \code{NULL} or the matrix of covariates one wishes to control
#' for, of size \code{nobs x ncovar}.
#' @param SNP_index \code{NULL} or the index vector of SNP to be tested. See
#' the 'Details' section.
#' @param alpha The significance level at which the FWER is controlled. Default
#' value is 0.05.
#' @param global.test Specifies wether the global null hypothesis should be
#' tested. Default value is \code{TRUE}. See the 'Details' section.
#' @param verbose Report information on progress. Default value is \code{TRUE}
#'
#' @details The testing is performed on the cluster tree given by \code{dendr}.
#' If the SNP data matrix was divided (e.g. by chromosome), and clustered
#' separately, the user must provide the argument \code{SNP_index}, to specify
#' which part of the data is being tested.
#'
#' Testing starts at the highest level, which includes all variables specified
#' by \code{SNP_index}, and moves down the cluster tree. It stops when a cluster's
#' null hypothesis cannot be rejected anymore. The smallest, still significant
#' clusters will be returned.
#'
#' By default the parameter \code{global.test = TRUE}, which means that first
#' the global null hypothesis is tested. If the data is divided (e.g. by
#' chromosome), and clustered separately, this parameter can be set to
#' \code{FALSE} once the global null has been rejected. This helps save time.
#'
#' @return A list of significant SNP groups with the following components:
#' \item{SNP_index}{The indeces of the SNPs in the group}
#' \item{pval}{The p-value of the SNP group}
#'
#' @examples
#' library(MASS)
#' x <- mvrnorm(60,mu = rep(0,60), Sigma = diag(60))
#' beta <- rep(0,60)
#' beta[c(5,9,3)] <- 1
#' y <- x %*% beta + rnorm(60)
#' dendr <- cluster.snp(x = x, method = "average")
#' res.multisplit <- multisplit(x, y)
#' sign.clusters <- test.hierarchy(x, y, dendr, res.multisplit)
#'
#' @references Buzdugan, L. et al. (2015), Assessing statistical significance in
#' predictive genome-wide association studies
#'
#' @name test.hierarchy
#' @aliases test.hierarchy

test.hierarchy <- function(x, y, dendr, res.multisplit, covar = NULL,
                           SNP_index = NULL, alpha = 0.05, global.test = TRUE,
                           verbose = TRUE) {
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
  # check if we need to cluster a subset only
  if (!is.null(SNP_index)) {
    if (length(SNP_index) > ncol(x)) {
      stop("The length of the SNP indeces is greater than the total number of
           SNPs")
    }
    if (length(which(SNP_index > ncol(x)) > 0)) {
      stop("There is no SNP with such index")
    }
    if (length(which(SNP_index < 1) > 0)) {
      stop("The SNP indeces must be positive integers")
    }
  }
  # compute the global p-value
  if (global.test == TRUE) {
      if (verbose) {
        message("Testing the global null hypothesis..")
      }
      ret.global <- comp.cluster.pval(seq(1,ncol(x)), res.multisplit, x, y,
                                      covar, is.global = TRUE)
  } else {
      ret.global <- list("pval" = 0)
  }
  # check if the global p-value is significant
  if (ret.global$pval <= alpha) {
    if (verbose) {
      message("The global null hypothesis was rejected")
    }
    pval.min <- ret.global$pval
    # check if we are testing a subset
    if (!is.null(SNP_index)) {
      if (verbose) {
        message("Testing a subset..")
      }
      ret.subset <- comp.cluster.pval(SNP_index, res.multisplit, x, y, covar)
      # redefine the maximum p-value
      if (ret.subset$pval > pval.min) {
        pval.min <- ret.subset$pval
      } else {
        ret.subset$pval <- pval.min
      }
      if (ret.subset$pval <= alpha) {
        if (verbose) {
          message("The null hypothesis of the subset was rejected")
          message("Testing the hierarchy of the subset..")
        }
        cluster.q <- vector('list',0)
        genotype.hclust.order <- order.dendrogram(dendr)
        cluster.q[[1]] <- list('clust' = dendr,
                               'SNP_index' = genotype.hclust.order,
                               'pval' = ret.subset$pval,
                               'pval.min' = pval.min)
        if (verbose == TRUE) {
          signif.cluster.list <- iterative.DFS(cluster.q, res.multisplit, x, y,
                                               alpha,SNP_index, covar)
        } else {
          signif.cluster.list <- suppressMessages(iterative.DFS(cluster.q,
                                                                res.multisplit,
                                                                x, y,
                                                                alpha, SNP_index,
                                                                covar))
        }
      } else {
        if (verbose) {
          message("The null hypothesis of the subset was not rejected")
          message("Testing stops.")
        }
        return(ret.global)
      }
    } else {
      if (verbose) {
        message("Testing the hierarchy")
      }
      cluster.q <- vector('list',0)
      genotype.hclust.order <- order.dendrogram(dendr)
      cluster.q[[1]] <- list('clust' = dendr,
                             'SNP_index' = genotype.hclust.order,
                             'pval' = ret.global$pval,
                             'pval.min' = pval.min)
      SNP_index <- seq(1,ncol(x))
      if (verbose == TRUE) {
        signif.cluster.list <- iterative.DFS(cluster.q, res.multisplit, x, y,
                                             alpha, seq(1,ncol(x)), covar)
      } else {
        signif.cluster.list <- suppressMessages(iterative.DFS(cluster.q,
                                                              res.multisplit,
                                                              x, y,
                                                              alpha,
                                                              seq(1,ncol(x)),
                                                              covar))
      }
    }
    return(signif.cluster.list)
  } else {
    if (verbose) {
      message("The global null hypothesis was accepted")
    }
    return(NULL)
  }
}
