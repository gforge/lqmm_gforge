
#' @rdname lqmm_predictions
predint.lqmm <- function(object, level = 0, alpha = 0.05, R = 50, seed = round(runif(1, 1, 10000)), newdata) {
  tau <- object$tau
  nq <- length(object$tau)
  p <- object$dim_theta[1]
  m <- object$dim_theta_z

  B <- boot(object, R = R, seed = seed)
  tmp <- object

  nobs <- object$nobs
  if (!missing(newdata)) {
    nobs <- nrow(newdata)
  }

  if (nq == 1) {
    yhat <- matrix(NA, nobs, R)
    for (i in 1:R) {
      tmp$theta <- B[i, 1:(p + m)]
      tmp$theta_x <- B[i, 1:p]
      tmp$theta_z <- B[i, (p + 1):(p + m)]
      tmp$scale <- B[i, (p + m + 1)]
      yhat[, i] <- predict(tmp, newdata = newdata, level = level)
    }
    LB <- apply(yhat, 1, quantile, probs = alpha / 2)
    UB <- apply(yhat, 1, quantile, probs = 1 - alpha / 2)
    ans <- data.frame(yhat = predict(object, newdata = newdata, level = level), lower = LB, upper = UB, SE = apply(yhat, 1, sd))
  } else {
    ans <- list()
    for (j in 1:nq) {
      tmp$tau <- tau[j]
      yhat <- matrix(NA, nobs, R)
      for (i in 1:R) {
        tmp$theta <- B[i, 1:(p + m), j]
        tmp$theta_x <- B[i, 1:p, j]
        tmp$theta_z <- B[i, (p + 1):(p + m), j]
        tmp$scale <- B[i, (p + m + 1), j]
        yhat[, i] <- predict(tmp, newdata = newdata, level = level)
      }
      LB <- apply(yhat, 1, quantile, probs = alpha / 2)
      UB <- apply(yhat, 1, quantile, probs = 1 - alpha / 2)
      ans[[j]] <- data.frame(yhat = predict(tmp, newdata = newdata, level = level), lower = LB, upper = UB, SE = apply(yhat, 1, sd))
    }
    names(ans) <- format(tau, digits = 4)
  }

  return(ans)
}
