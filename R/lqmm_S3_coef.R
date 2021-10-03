#' Extract LQMM Coefficients
#'
#' \code{coef} extracts model coefficients from \code{lqmm} objects.
#'
#'
#' @param object a fitted object of \code{\link{class}} "lqmm".
#' @param \dots not used.
#' @return a vector for single quantiles or a matrix for multiple quantiles.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}} \code{\link{summary.lqmm}}
#' @keywords coefficients
coef.lqmm <- function(object, ...) {
  tau <- object$tau
  nq <- length(tau)
  ans <- object$theta_x

  if (nq == 1) {
    names(ans) <- object$nn
  }

  return(ans)
}
