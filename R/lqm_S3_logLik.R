#' Extract Log-Likelihood
#' 
#' \code{logLik.lqm} extracts the log-likelihood of a fitted LQM.
#' 
#' 
#' @param object an object of \code{\link{class}} "lqm".
#' @param \dots not used.
#' @author Marco Geraci
#' @seealso \code{\link{lqm}} \code{\link{AIC}}
#' @keywords models
logLik.lqm <- function(object, ...) {
  tdf <- object$edf + 1
  tau <- object$tau
  nq <- length(tau)

  if (nq == 1) {
    ans <- object$logLik
  } else {
    ans <- NULL
    for (i in 1:nq) ans <- c(ans, object[[i]]$logLik)
    names(ans) <- as.character(format(tau, digits = 4))
  }

  attr(ans, "nobs") <- object$nobs
  attr(ans, "df") <- tdf
  attr(ans, "class") <- "logLik"

  return(ans)
}
