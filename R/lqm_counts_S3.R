coef.lqm.counts <- function(object, ...) {
  tau <- object$tau
  nq <- length(tau)
  ans <- object$theta

  if (nq == 1) {
    names(ans) <- object$term.labels
  }

  return(ans)
}

predict.lqm.counts <- function(object, newdata, na.action = na.pass, ...) {
  tau <- object$tau

  if (missing(newdata)) {
    yhat <- drop(object$x %*% object$theta)
  } else {
    objt <- terms(object)
    Terms <- delete.response(objt)
    m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    yhat <- drop(x %*% object$theta)
  }

  return(yhat)
}

residuals.lqm.counts <- function(object, ...) {
  ans <- as.numeric(object$y) - predict(object)
  return(ans)
}

print.lqm.counts <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  tau <- x$tau
  nq <- length(tau)
  cat("Call: ")
  dput(x$call)
  cat("\n")
  if (nq == 1) {
    cat(paste("Quantile", tau, "\n"))
    cat("\n")
    cat("Fixed effects:\n")
    printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
  } else {
    NULL
  }
}
