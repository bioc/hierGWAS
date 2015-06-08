# This function computes the adjusted aggregated p value of a cluster
# Input: pval.vector - the vector of p values for a cluster
# Output: adjusted.pval - the adjusted aggregated p value of a cluster

adj.pval <- function (ret.vector, B) {
  # define the sequence of gamma values
  gamma.min <- 0.05
  gamma.step <- 0.01
  gamma.seq <- seq(gamma.min,1,gamma.step)

  pval.vector <- vector("numeric",B)
  for (i in 1:B) {
    pval.vector[i] <- ret.vector[[i]]$pval
  }

  # compute the empirical quantile vector
  quantile.vector <- vector("numeric", length = length(gamma.seq))
  for (g in 1:length(gamma.seq)) {
    quantile.vector[g] <- min(1, quantile(pval.vector/gamma.seq[g], gamma.seq[g]
                                          , na.rm = TRUE))
  }

  # compute the adjusted p value
  adjusted.pval <- min(1, (1 - log(gamma.min))*min(quantile.vector))
  return(list("pval" = adjusted.pval))
}
