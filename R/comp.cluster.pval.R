# This function computes the adjusted p value of a cluster of any size
# Inputs: SNP_index - the variables that belong to the cluster
#         multisplit.out - the output of the multisplit function
#         x - the predictor matrix
#         y - the response
#         covar - matrix of covariates
# Output: adjusted.pval - the adjusted p value of a cluster


comp.cluster.pval <- function (SNP_index, multisplit.out, x, y, covar = NULL,
                               is.global = FALSE) {
  # get the samples and coefficients
  out.sample.matrix <- multisplit.out$out.sample
  lasso.coeff.matrix <- multisplit.out$sel.coeff
  B <- nrow(lasso.coeff.matrix)
  ret.vector <- vector("list", length = B)
  s <- ncol(lasso.coeff.matrix)
  # compute the p value for each split
  for (b in 1:B) {
    scl <- lasso.coeff.matrix[b,]
    scl <- scl[!is.na(scl)]
    common.SNP_index <- intersect(SNP_index,scl)
    common.SNP_index.index <- match(common.SNP_index,scl)
    if (length(common.SNP_index) == 0) {
      ret.vector[[b]]$pval <- 1
    } else {
      c <- length(common.SNP_index)
      cf <- s/c
      #show(cf)
      x.b <- x[out.sample.matrix[b, ], scl]
      y.b <- y[ out.sample.matrix[b, ]]
      if (is.null(covar)) {
        ret.vector[[b]] <- test.snp(x.b, y.b,
                                    cluster.index = common.SNP_index.index,
                                    is.global = is.global)
      } else {
        covar <- data.matrix(covar)
        ret.vector[[b]] <- test.snp(x.b, y.b,
                                    cluster.index = common.SNP_index.index,
                                    is.global = is.global,
                                    covar = covar[out.sample.matrix[b, ],])
      }

      pval <- ret.vector[[b]]$pval
      pval <- min(pval*cf,1)
      ret.vector[[b]]$pval <- pval
    }
  }
  # compute the adjusted p value
  adjusted.ret <- adj.pval(ret.vector, B)
  return(adjusted.ret)
}


