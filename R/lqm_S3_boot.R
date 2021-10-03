boot.lqm <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE) {
  set.seed(seed)
  tau <- object$tau
  nq <- length(tau)
  obsS <- replicate(R, sample(1:object$nobs, replace = TRUE))
  npars <- object$dim_theta

  control <- object$control
  control$verbose <- FALSE

  if (nq == 1) {
    bootmat <- matrix(NA, R, npars)
    colnames(bootmat) <- object$term.labels
    FIT_ARGS <- list(theta = object$theta, tau = tau, control = control)
    for (i in 1:R) {
      a <- table(obsS[, i])
      s <- as.numeric(names(a))
      FIT_ARGS$x <- as.matrix(object$x[s, ])
      FIT_ARGS$y <- object$y[s]
      FIT_ARGS$weights <- as.numeric(a)
      FIT_ARGS$theta <- if (!startQR) lm.wfit(x = FIT_ARGS$x, y = FIT_ARGS$y, w = FIT_ARGS$weights)$coefficients
      fit <- try(do.call(lqm.fit.gs, FIT_ARGS), silent = TRUE)
      if (!inherits(fit, "try-error")) bootmat[i, ] <- fit$theta
    }
  } else {
    bootmat <- array(NA, dim = c(R, npars, nq), dimnames = list(NULL, object$term.labels, paste("tau = ", format(tau, digits = 4), sep = "")))
    FIT_ARGS <- list(control = control)
    for (i in 1:R) {
      a <- table(obsS[, i])
      s <- as.numeric(names(a))
      for (j in 1:nq) {
        FIT_ARGS$x <- as.matrix(object$x[s, ])
        FIT_ARGS$y <- object$y[s]
        FIT_ARGS$weights <- as.numeric(a)
        FIT_ARGS$theta <- if (startQR) object[[j]]$theta else lm.wfit(x = FIT_ARGS$x, y = FIT_ARGS$y, w = FIT_ARGS$weights)$coefficients
        FIT_ARGS$tau <- tau[j]
        fit <- try(do.call(lqm.fit.gs, FIT_ARGS), silent = TRUE)
        if (!inherits(fit, "try-error")) bootmat[i, , j] <- fit$theta
      }
    }
  }

  class(bootmat) <- "boot.lqm"
  attr(bootmat, "tau") <- tau
  attr(bootmat, "estimated") <- object$theta
  attr(bootmat, "R") <- R
  attr(bootmat, "seed") <- seed
  attr(bootmat, "npars") <- npars
  attr(bootmat, "indices") <- obsS
  attr(bootmat, "rdf") <- object$rdf

  return(bootmat)
}



#' Summary for a \code{boot.lqm} Object
#' 
#' Summary method for class \code{boot.lqm}.
#' 
#' 
#' @param object an object of \code{\link{class}} \code{lqm}.
#' @param alpha numeric value for the interval confidence level
#' (\code{1-alpha}).
#' @param digits a non-null value for digits specifies the minimum number of
#' significant digits to be printed in values.
#' @param \dots not used.
#' @author Marco Geraci
#' @seealso \code{\link{boot.lqm}}, \code{\link{lqm}},
#' @keywords summary bootstrap
summary.boot.lqm <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), ...) {
  tau <- attr(object, "tau")
  nq <- length(tau)

  est <- attr(object, "estimated")
  npars <- attr(object, "npars")
  rdf <- attr(object, "rdf")
  R <- attr(object, "R")


  nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound", "Pr(>|t|)")

  if (nq == 1) {
    bias <- est - apply(as.matrix(object), 2, mean)
    Cov <- cov(as.matrix(object))
    stds <- sqrt(diag(Cov))
    lower <- est + qt(alpha / 2, R - 1) * stds
    upper <- est + qt(1 - alpha / 2, R - 1) * stds
    tP <- 2 * pt(-abs(est / stds), R - 1)
    ans <- cbind(est, bias, stds, lower, upper, tP)
    colnames(ans) <- nn
    printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
  } else {
    bias <- est - apply(object, 3, colMeans)
    Cov <- apply(object, 3, function(x) cov(as.matrix(x)))
    if (npars == 1) Cov <- matrix(Cov, nrow = 1)
    stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
    lower <- est + qt(alpha / 2, R - 1) * stds
    upper <- est + qt(1 - alpha / 2, R - 1) * stds
    tP <- 2 * pt(-abs(est / stds), R - 1)
    for (i in 1:nq) {
      if (npars == 1) {
        ans <- c(est[i], bias[i], stds[i], lower[i], upper[i], tP[i])
        ans <- matrix(ans, nrow = 1)
      } else {
        ans <- cbind(est[, i], bias[, i], stds[, i], lower[, i], upper[, i], tP[, i])
      }
      rownames(ans) <- rownames(est)
      colnames(ans) <- nn
      cat(paste("tau = ", tau[i], "\n", sep = ""))
      printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
      cat("\n")
    }
  }
}
