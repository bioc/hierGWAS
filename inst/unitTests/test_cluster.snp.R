test_cluster.snp <- function() {
  library(MASS)
  library(RUnit)
  d <- mvrnorm(50,mu = rep(0,50), Sigma = diag(50))
  x <- mvrnorm(50,mu = rep(0,100), Sigma = diag(100))
  SNP_index <- seq(1,10)

  dendr <- cluster.snp(x)
  checkEquals(length(order.dendrogram(dendr)),100)

  dendr <- cluster.snp(d = d, method = "complete")
  checkEquals(length(order.dendrogram(dendr)),50)

  dendr <- cluster.snp(x = x, d = d, method = "complete")
  checkEquals(length(order.dendrogram(dendr)),50)

  dendr <- cluster.snp(x, method = "complete", SNP_index = SNP_index)
  checkEquals(length(order.dendrogram(dendr)),10)

  dendr <- cluster.snp(d = d, method = "complete", SNP_index = SNP_index)
  checkEquals(length(order.dendrogram(dendr)),50)

  dendr <- cluster.snp(x = x, d = d, method = "complete", SNP_index = SNP_index)
  checkEquals(length(order.dendrogram(dendr)),50)
}
