# MEL: maximum estimated likelihood estimator
# Copyright (C) Peter J. Rousseeuw, Andreas Christmann, 2001
# Please note, that there is NO WARRANTY, 18/JUN/2001

"MEL" <- function(x, y, delta=0.01, epsilon=1.0E-6, maxit=100) {
  pihat <- max(delta, min(1-delta, mean(y)))
  delta0 <- (pihat*delta) / (1+delta)
  delta1 <- (1+pihat*delta) / (1+delta)
  ytilde <- delta0 * (1-y) + delta1 * y
  response <- cbind(ytilde, 1-ytilde)
  suppressWarnings(outMEL <- glm(response ~ x, family=binomial, 
          control=glm.control(epsilon=epsilon, maxit=maxit)))
  beta <- coef(outMEL)
  ### construct output object
  output <- list(MEL=beta, outMEL=outMEL)
  class(output) <- "MEL"
  output$call <- match.call()
  output
}



