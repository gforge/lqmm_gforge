#' Residuals from an LQM Objects
#' 
#' This function computes the residuals from a fitted linear quantile model.
#' 
#' 
#' @aliases residuals.lqm residuals.lqm.counts
#' @param object an \code{lqm} or \code{lqm.counts} object.
#' @param \dots not used.
#' @return a vector or matrix of residuals.
#' @author Marco Geraci
#' @seealso \code{\link{lqm}}, \code{\link{lqm.counts}},
#' \code{\link{predict.lqm}}, \code{\link{coef.lqm}}
#' @keywords residuals
residuals.lqm <- function(object, ...) {
  ans <- as.numeric(object$y) - predict(object)
  return(ans)
}
