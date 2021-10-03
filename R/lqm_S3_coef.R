coef.lqm <- function(object, ...) {
  tau <- object$tau
  nq <- length(tau)
  ans <- object$theta

  if (nq == 1) {
    names(ans) <- object$term.labels
  }

  return(ans)
}
