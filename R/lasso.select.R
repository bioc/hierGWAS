# This function runs the lasso on a random data split
# Input: x - the matrix of predictors
#        y - the vector of responses
#        ncovae - the number of covariates
#        pen - the penalty (covariates are not penalized)
#        family - continuous or binary response
#        i - the iteration number
# Output: out.sample - the second half sample
#         selected.lasso.coeff - the indeces of the selected lasso coefficients

lasso.select <- function(x, y, ncovar, pen, family, i) {
  # get the in.sample and out.sample
  no.samples <- length(y)
  no.snps <- dim(x)[2] - ncovar
  in.sample.length <- floor(no.samples/2)
  in.sample <- sort(sample(seq(1:no.samples),in.sample.length))
  out.sample <- sort(setdiff(seq(1:no.samples),in.sample))

  # perform lasso on the in.sample
  lasso.output <- glmnet( x[in.sample,], y[in.sample],
                          family = family, standardize = FALSE,
                          intercept = TRUE, penalty.factor = pen)

  # save the first n/6 regression coefficients that enter the lasso path
  beta.matrix <- lasso.output$beta
  lasso.coeff.path <- beta.matrix@i + 1
  unique.lasso.coeff <- unique(lasso.coeff.path)
  # check if we have covariates
  if (ncovar > 0) {
    unique.lasso.coeff <- setdiff(unique.lasso.coeff,seq(1,ncovar))
    selected.lasso.coeff <- unique.lasso.coeff[1 : (no.samples/6)] - ncovar
  } else {
    selected.lasso.coeff <- unique.lasso.coeff[1 : (no.samples/6)]
  }
  return(list("out.sample" = out.sample, "sel.coeff" = selected.lasso.coeff))
}
