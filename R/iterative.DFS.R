iterative.DFS <- function(cluster.queue, multisplit.out, x, y,
                          alpha, subset.index, covar) {
  # the list that stores the clusters that were found significant
  signif.cluster.list <- vector('list',0)
  while (length(cluster.queue) > 0) {
    # the total number of elements in the cluster
    clust <- cluster.queue[[1]]$clust
    member.SNP_index <- cluster.queue[[1]]$SNP_index
    adjusted.pval <- cluster.queue[[1]]$pval
    pval.min <- cluster.queue[[1]]$pval.min
    if (adjusted.pval < pval.min) {
      adjusted.pval <- pval.min
    } else {
      pval.min <- adjusted.pval
    }
    total.el <- length(member.SNP_index)
    h <- attr(clust,'height')
    # check if the cluster is a leaf
    if (total.el <= 1) {
      # put the significant leaf in the result list
      cluster.queue <- cluster.queue[-1]
      signif.cluster.list[[length(signif.cluster.list)+1]] <- list(
        'label' = subset.index[member.SNP_index], 'pval' = adjusted.pval)
    }
    else {
      # get the subclusters of the current cluster
      subclust = cut(clust, h = attr(clust,'height'))$lower
      # the number of elements in the 1st subcluster
      no.el1 <- attr(subclust[[1]],'members')
      # the members that belong to the 2 subclusters
      SNP_index.child1 <- member.SNP_index[1:no.el1]
      SNP_index.child2 <- member.SNP_index[(no.el1+1):total.el]
      no.el1 <- length(SNP_index.child1)
      no.el2 <- length(SNP_index.child2)
      # check the children of this cluster
      mess <- paste("Tested a group with",no.el1 + no.el2,"SNPs.")
      message(mess)
      message("The group is significant.")
      mess <- paste("Testing the children of the group:")
      message(mess)
      mess <- paste("Child 1 with",no.el1,"SNPs and child 2 with",no.el2)
      message(mess)
      adjusted.ret1 <- comp.cluster.pval(subset.index[SNP_index.child1],
                                         multisplit.out, x, y, covar)
      adjusted.ret2 <- comp.cluster.pval(subset.index[SNP_index.child2],
                                         multisplit.out, x, y, covar)
      adjusted.pval1 <- adjusted.ret1$pval
      adjusted.pval2 <- adjusted.ret2$pval
      flag.child1 = (adjusted.pval1 <= alpha)
      flag.child2 = (adjusted.pval2 <= alpha)
      # the number of elements that are in the cluster queue
      ll <- length(cluster.queue)
      # neither of the children are significant
      # put the current cluster in the result list and stop
      if ((flag.child1 == 0) && (flag.child2 == 0)) {
        cluster.queue <- cluster.queue[-1]
        signif.cluster.list[[length(signif.cluster.list)+1]] <- list(
          'label' = subset.index[member.SNP_index], 'pval' = adjusted.pval)
      }
      # child 1 is significant
      # remove the current cluster from the queue, add child 1 to the queue
      if ((flag.child1 == 1) && (flag.child2 == 0)) {
        cluster.queue <- cluster.queue[-1]
        cluster.queue[[ll]] <- list('clust' = subclust[[1]],
                                    'SNP_index' = SNP_index.child1,
                                    'pval' = adjusted.pval1,
                                    'pval.min' = pval.min)
      }
      # child 2 is significant
      # remove the current cluster from the queue, add child 2 to the queue
      if ((flag.child1 == 0) && (flag.child2 == 1)) {
        cluster.queue <- cluster.queue[-1]
        cluster.queue[[ll]] <- list('clust' = subclust[[2]],
                                    'SNP_index' = SNP_index.child2,
                                    'pval' = adjusted.pval2,
                                    'pval.min' = pval.min)
      }
      # both children are significant
      # remove the current cluster from the queue, add both children to it
      if ((flag.child1 == 1) && (flag.child2 == 1)) {
        cluster.queue <- cluster.queue[-1]
        cluster.queue[[ll]] <- list('clust' = subclust[[1]],
                                    'SNP_index' = SNP_index.child1,
                                    'pval' = adjusted.pval1,
                                    'pval.min' = pval.min)
        cluster.queue[[ll+1]] <- list('clust' = subclust[[2]],
                                      'SNP_index' = SNP_index.child2,
                                      'pval' = adjusted.pval2,
                                      'pval.min' = pval.min)
      }
    }
  }
  # return the list of significant clusters
  return(signif.cluster.list)
}

