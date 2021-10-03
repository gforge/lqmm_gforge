#' Extract LQM Coefficients
#' 
#' \code{coef} extracts model coefficients from \code{lqm}, \code{lqm.counts}
#' objects.
#' 
#' 
#' @aliases coef.lqm coef.lqm.counts
#' @param object an \code{lqm} or \code{lqm.counts} object.
#' @param \dots not used.
#' @return a vector for single quantiles or a matrix for multiple quantiles.
#' @author Marco Geraci
#' @seealso \code{\link{lqm}} \code{\link{summary.lqm}}
#' \code{\link{lqm.counts}}
#' @keywords coefficients
coef.lqm <- function(object, ...) {
  tau <- object$tau
  nq <- length(tau)
  ans <- object$theta

  if (nq == 1) {
    names(ans) <- object$term.labels
  }

  return(ans)
}
