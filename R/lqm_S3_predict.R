#' Predictions from LQM Objects
#'
#' This function computes predictions based on fitted linear quantile model.
#'
#' @param object an lqm or lqm.counts object.
#' @param newdata an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param interval logical flag. If TRUE, bootstrap percentile intervals for predictions are provided. This argument is for lqm objects only.
#' @param level confidence level. This argument is for lqm objects only.
#' @param na.action function determining what should be done with missing values in newdata. The default is to predict NA.
#' @param ... further arguments passed to [boot.lqm].
#'
#' @return a vector or a matrix or an array of predictions.
#' @export
#' @author Marco Geraci
#' @seealso \code{\link{residuals.lqm}}, \code{\link{residuals.lqm.counts}},
#' \code{\link{lqm}}, \code{\link{lqm.counts}}, \code{\link{coef.lqm}},
#' \code{\link{boot.lqm}}
#' @keywords prediction
#' @aliases predict.lqm predict.lqm.counts
predict.lqm <- function(object, newdata, interval = FALSE,
                        level = 0.95, na.action = na.pass, ...) {
  tau <- object$tau
  nq <- length(tau)
  lp <- (1 - level) / 2
  up <- 1 - lp

  qrow <- function(x, p) apply(x, 1, function(y) quantile(y, p))

  if (missing(newdata)) {
    X <- object$x
  } else {
    objt <- terms(object)
    Terms <- delete.response(objt)
    m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  }

  if (nq == 1) {
    yhat <- drop(X %*% object$theta)
    if (interval) {
      lin_pred <- X %*% t(boot.lqm(object, ...))
      low <- qrow(lin_pred, lp)
      upp <- qrow(lin_pred, up)
      yhat <- cbind(yhat, low, upp)
      colnames(yhat) <- c("fitted", "lower", "higher")
    }
  } else {
    yhat <- matrix(X %*% as.matrix(object$theta), ncol = nq)
    if (interval) {
      B <- boot.lqm(object, ...)
      Rboot <- attr(B, "R")
      lin_pred <- matrix(apply(B, 3, function(b, x) x %*% as.matrix(t(b)), x = X), ncol = nq)
      low <- matrix(apply(lin_pred, 2, function(x, R, p) qrow(matrix(x, ncol = R), p), R = Rboot, p = lp), ncol = nq)
      upp <- matrix(apply(lin_pred, 2, function(x, R, p) qrow(matrix(x, ncol = R), p), R = Rboot, p = up), ncol = nq)
      res <- array(NA, dim = c(nrow(X), 3, nq), dimnames = list(NULL, c("fitted", "lower", "higher"), format(tau, 4)))
      for (i in 1:nq) res[, , i] <- cbind(yhat[, i], low[, i], upp[, i])
      yhat <- res
    }
  }

  return(yhat)
}
