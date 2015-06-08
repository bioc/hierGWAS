# return the R2 for a linear or logistic regression model

return.r2 <- function (x, y, covar = NULL) {

  # check if the response is binary or continuous
  if ((min(y) == 0) & (max(y) == 1) & (length(unique(y)) == 2)) {
    is.binary = TRUE
  } else {
    is.binary = FALSE
  }
  if (is.binary) {
    if (is.null(covar)) {
      design.mat <- x
    } else {
      design.mat <- data.matrix(cbind(covar, x))
    }
    r2 <- NagelkerkeR2(glm(y~design.mat, family = binomial(link="logit")))$R2
  } else {
    if (is.null(covar)) {
      design.mat <- x
    } else {
      design.mat <- data.matrix(cbind(covar, x))
    }
    r2 <- summary(lm(y~design.mat))$adj.r.squared
  }
  return(r2)
}
