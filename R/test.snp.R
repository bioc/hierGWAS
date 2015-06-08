# This function performs a full or partial F test, or full or partial LRT
# Input: x - the matrix of predictors
#        y - the vector of responses
#        cluster.index - the set of variables to be tested
#        is.global - the flag that specifies if it is the global cluster
# Output: the p value of the test

test.snp <- function (x, y, cluster.index, is.global, covar = NULL) {

  # check if the response is binary or continuous
  if ((min(y) == 0) & (max(y) == 1) & (length(unique(y)) == 2)) {
    is.binary = TRUE
  } else {
    is.binary = FALSE
  }
  if (is.binary) {
    if (is.global){
      if (is.null(covar)) {
        design.mat <- x
        control.mat <- rep(1,length(y))
      } else {
        design.mat <- data.matrix(cbind(covar, x))
        control.mat <- data.matrix(covar)
      }
      glm.1 <- MEL(control.mat, y, maxit = 100)$outMEL
      glm.2 <- MEL(design.mat, y, maxit = 100)$outMEL
      a <- anova(glm.1,glm.2, test = "Chisq")
      pval <- a$"Pr(>Chi)"[2]
    }
    else{
      compl.cluster.index <- setdiff( (1 : ncol(x)), cluster.index)
      if (is.null(covar)) {
        design.mat <- x
        design.mat.partial <- x[,compl.cluster.index]
      } else {
        design.mat <- data.matrix(cbind(covar, x))
        design.mat.partial <- data.matrix(cbind(covar, x[,compl.cluster.index]))
      }
      glm.1 <- MEL(design.mat.partial, y, maxit = 100)$outMEL
      glm.2 <- MEL(design.mat, y, maxit = 100)$outMEL
      a <- anova(glm.1, glm.2, test = "Chisq")
      pval <- a$"Pr(>Chi)"[2]
    }
  } else {
    if (is.global){
      if (is.null(covar)) {
        pval <- anova(lm(y~x))$"Pr(>F)"[1]
      } else {
        design.mat <- data.matrix(cbind(covar, x))
        control.mat <- data.matrix(covar)
        pval <- anova(lm(y~control.mat),lm(y~design.mat), test = "F")$P[2]
      }
    }
    else{
      compl.cluster.index <- setdiff( (1 : ncol(x)), cluster.index)
      if (is.null(covar)) {
        design.mat <- x
        design.mat.partial <- x[,compl.cluster.index]
      } else {
        design.mat <- data.matrix(cbind(covar, x))
        design.mat.partial <- data.matrix(cbind(covar, x[,compl.cluster.index]))
      }
      # do a partial F test
      pval <- anova(lm(y~design.mat.partial),
                    lm(y~design.mat), test = "F")$P[2]
    }
  }

  return(list("pval" = pval))
}

