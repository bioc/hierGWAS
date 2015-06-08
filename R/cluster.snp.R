#' Hierarchical Clustering of SNP Data
#'
#' Clusters SNPs hierachically.
#'
#' @param x The SNP data matrix of size \code{nobs x nvar}. Default value is
#' \code{NULL}
#' @param d \code{NULL} or a dissimilarity matrix. See the 'Details' section.
#' @param method The agglomeration method to be used. This should be
#' (an unambiguous abbreviation of) one of \code{"ward.D"}, \code{"ward.D2"},
#' \code{"single"}, \code{"complete"}, \code{"average"} (= UPGMA),
#' \code{"mcquitty"} (= WPGMA), \code{"median"} (= WPGMC) or \code{"centroid"}
#' (= UPGMC). See \code{\link[fastcluster]{hclust}} for details.
#' @param SNP_index \code{NULL} or the index vector of SNPs to be clustered. See
#' the 'Details' section.
#'
#' @details The SNPs are clustered using \code{\link[fastcluster]{hclust}},
#' which performs a hierarchical cluster analysis using a set of dissimilarities
#' for the \code{nvar} objects being clustered. There are 3 possible scenarios.
#'
#' If \code{d = NULL}, \code{x} is used to compute the dissimilarity matrix.
#' The dissimilarity measure between two SNPs is \code{1 - LD} (Linkage
#' Disequilibrium), where \code{LD} is defined as the square of the Pearson
#' correlation coefficient. If \code{SNP_index = NULL}, all \code{nvar} SNPs will
#' be clustered; otherwise only the SNPs with indices specified by \code{SNP_index}
#' will be considered.
#'
#' If the user wishes to use a different dissimilarity measure, \code{d} needs
#' to be provided. \code{d} must be either a square matrix of size
#' \code{nvar x nvar}, or an object of class \code{\link{dist}}. If \code{d} is
#' provided, \code{x} and \code{SNP_index} will be ignored.
#'
#' @return An object of class \code{\link{dendrogram}} which describes the tree
#' produced by the clustering algorithm \code{\link[fastcluster]{hclust}}.
#'
#' @examples
#' library(MASS)
#' x <- mvrnorm(60,mu = rep(0,60), Sigma = diag(60))
#' clust.1 <- cluster.snp(x = x, method = "average")
#' SNP_index <- seq(1,10)
#' clust.2 <- cluster.snp(x = x, method = "average", SNP_index = SNP_index)
#' d <- dist(x)
#' clust.3 <- cluster.snp(d = d, method = "single")
#' @name cluster.snp


cluster.snp <- function(x = NULL, d = NULL, method = "average",
                        SNP_index = NULL) {

  if (!is.null(x)) {
    if (!is.numeric(x)) {
      stop("'x' must be numeric")
    }
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
  }
  # check if we need to cluster a subset only
  if (!is.null(SNP_index) && !is.null(x)) {
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
    x <- x[,SNP_index]
  }
  # compute the distance matrix, if it was not provided
  if (is.null(d)) {
    d <- 1 - abs(cor(x))^2
  } else {
    if (!is.numeric(d)) {
      stop("'d' must be numeric")
    }
  }
  dist.matrix <- as.dist(d)
  # cluster the SNPs
  x.hclust <- hclust(dist.matrix, method = method)
  x.dendr <- as.dendrogram(x.hclust)
  return(x.dendr)
}
