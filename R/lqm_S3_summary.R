#' Summary for an \code{lqm} Object
#' 
#' Summary method for class \code{lqm}.
#' 
#' \code{print.summary.lqm} formats the coefficients, standard errors, etc. and
#' additionally gives `significance stars'.
#' 
#' @param object an object of \code{\link{class}} \code{lqm}
#' @param method specifies the method used to compute standard errors: "boot"
#' for bootstrap (default), "nid" for large sample approximations under
#' \emph{nid} assumptions.
#' @param alpha significance level.
#' @param covariance logical flag. If \code{TRUE} the covariance matrix of the
#' bootstrap estimates is provided.
#' @param \dots see \code{\link{boot.lqm}} for additional arguments.
#' @return an object of class \code{summary.lqm}. The function
#' \code{summary.lqm} computes and returns a list of summary statistics of the
#' fitted linear quantile mixed model given in \code{object}, using the
#' components (list elements) from its argument, plus
#' 
#' \item{Cov}{the covariance matrix obtained from the bootstrapped estimates
#' (if \code{covariance = TRUE}).} \item{tTable}{a matrix with estimates,
#' standard errors, etc.}
#' @author Marco Geraci
#' @seealso \code{\link{print.summary.lqm}} \code{\link{lqm}}
#' @source The code for the "nid" method has been adapted from the function
#' \code{summary.rq} in package \code{quantreg}. It depends on the function
#' \code{bandwidth.rq}.
#' 
#' Roger Koenker (2016). quantreg: Quantile Regression. R package version 5.29.
#' \url{https://CRAN.R-project.org/package=quantreg}
#' @keywords bootstrap standard errors
#' @examples
#' 
#' 
#' set.seed(12356)
#' n <- 200
#' p <- 1:3/4
#' test <- data.frame(x = runif(n,0,1))
#' test$y <- 30 + test$x + rnorm(n)
#' fit.lqm <- lqm(y ~ x, data = test, tau = p)
#' summary(fit.lqm, R = 50)
#' 
#' 
summary.lqm <- function(object, method = "boot", alpha = 0.05, covariance = FALSE, ...) {
  tau <- object$tau
  nq <- length(tau)
  theta <- object$theta
  x <- object$x
  y <- object$y
  n <- nrow(x)
  npars <- object$dim_theta
  rdf <- object$rdf

  nn <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")

  if (method == "boot") {
    B <- boot.lqm(object, ...)
    R <- attr(B, "R")
    if (nq == 1) {
      Cov <- cov(as.matrix(B))
      stds <- sqrt(diag(Cov))
      tP <- 2 * pt(-abs(theta / stds), R - 1)
      lower <- theta + qt(alpha / 2, R - 1) * stds
      upper <- theta + qt(1 - alpha / 2, R - 1) * stds
      ans <- cbind(theta, stds, lower, upper, tP)
      colnames(ans) <- nn
    } else {
      Cov <- apply(B, 3, function(x) cov(as.matrix(x)))
      if (npars == 1) Cov <- matrix(Cov, nrow = 1)
      stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
      tP <- 2 * pt(-abs(theta / stds), R - 1)
      lower <- theta + qt(alpha / 2, R - 1) * stds
      upper <- theta + qt(1 - alpha / 2, R - 1) * stds
      ans <- vector("list", nq)
      Cov.array <- array(NA, dim = c(npars, npars, nq))
      for (i in 1:nq) {
        if (npars == 1) {
          ans[[i]] <- matrix(c(theta[i], stds[i], lower[i], upper[i], tP[i]), nrow = 1)
          rownames(ans[[i]]) <- rownames(theta)
          colnames(ans[[i]]) <- nn
        } else {
          ans[[i]] <- cbind(theta[, i], stds[, i], lower[, i], upper[, i], tP[, i])
          rownames(ans[[i]]) <- rownames(theta)
          colnames(ans[[i]]) <- nn
        }
        Cov.array[, , i] <- matrix(Cov[, i], npars, npars)
      }
      Cov <- Cov.array
      dimnames(Cov) <- list(rownames(theta), rownames(theta), format(tau, digits = 4))
    }
  }
  if (method == "nid") {
    eps <- .Machine$double.eps^(2 / 3)
    h <- bandwidth.rq(tau, n, hs = TRUE)
    if (nq == 1) {
      if (tau + h > 1 || tau - h < 0) stop("bandwith is too large or tau is too close to 0 or 1")
      theta.up <- update(object, tau = tau + h, evaluate = TRUE)$theta
      # theta.up <- lqm.fit.gs(theta, x, y, object$weights, tau + h, object$control)$theta
      theta.low <- update(object, tau = tau - h, evaluate = TRUE)$theta
      # theta.low <- lqm.fit.gs(theta, x, y, object$weights, tau - h, object$control)$theta
      dyhat <- x %*% matrix(theta.up - theta.low, ncol = 1)
      dens <- pmax(0, (2 * h) / (dyhat - eps))
      fxxinv <- diag(npars)
      fxxinv <- backsolve(qr(sqrt(dens) * x)$qr[1:npars, 1:npars, drop = FALSE], fxxinv)
      fxxinv <- fxxinv %*% t(fxxinv)
      Cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
      dimnames(Cov) <- list(colnames(x), colnames(x))
      stds <- sqrt(diag(Cov))
      tP <- 2 * pt(-abs(theta / stds), object$rdf)
      lower <- theta + qt(alpha / 2, object$rdf) * stds
      upper <- theta + qt(1 - alpha / 2, object$rdf) * stds
      ans <- cbind(theta, stds, lower, upper, tP)
      colnames(ans) <- nn
    } else {
      if (any(tau + h > 1) || any(tau - h < 0)) stop("bandwith is too large or tau is too close to 0 or 1")
      theta.up <- update(object, tau = tau + h, evaluate = TRUE)$theta
      theta.low <- update(object, tau = tau - h, evaluate = TRUE)$theta
      dyhat <- x %*% (theta.up - theta.low)
      Cov.array <- array(NA, dim = c(npars, npars, nq))
      for (i in 1:nq) {
        dens <- pmax(0, (2 * h[i]) / (dyhat[, i] - eps))
        fxxinv <- diag(npars)
        fxxinv <- backsolve(qr(sqrt(dens) * x)$qr[1:npars, 1:npars, drop = FALSE], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        Cov <- tau[i] * (1 - tau[i]) * fxxinv %*% crossprod(x) %*% fxxinv
        Cov.array[, , i] <- Cov
      }
      stds <- apply(Cov.array, 3, function(x) sqrt(diag(x)))
      tP <- 2 * pt(-abs(theta / stds), object$rdf)
      lower <- theta + qt(alpha / 2, object$rdf) * stds
      upper <- theta + qt(1 - alpha / 2, object$rdf) * stds
      ans <- vector("list", nq)
      for (i in 1:nq) {
        if (npars == 1) {
          ans[[i]] <- matrix(c(theta[i], stds[i], lower[i], upper[i], tP[i]), nrow = 1)
          rownames(ans[[i]]) <- rownames(theta)
          colnames(ans[[i]]) <- nn
        } else {
          ans[[i]] <- cbind(theta[, i], stds[, i], lower[, i], upper[, i], tP[, i])
          rownames(ans[[i]]) <- rownames(theta)
          colnames(ans[[i]]) <- nn
        }
      }
      Cov <- Cov.array
      dimnames(Cov) <- list(rownames(theta), rownames(theta), format(tau, digits = 4))
    }
  }

  if (covariance) object$Cov <- Cov
  object$tTable <- ans

  class(object) <- "summary.lqm"
  return(object)
}



#' Print an \code{lqm} Summary Object
#' 
#' Print summary of an \code{lqm} object.
#' 
#' 
#' @param x a \code{summary.lqm} object.
#' @param \dots not used.
#' @author Marco Geraci
#' @seealso \code{\link{lqm}}, \code{\link{summary.lqm}}
#' @keywords print summary
print.summary.lqm <- function(x, ...) {
  tau <- x$tau
  nq <- length(tau)

  cat("Call: ")
  dput(x$call)

  if (nq == 1) {
    cat(paste("Quantile", tau, "\n"))
    printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
  } else {
    for (i in 1:nq) {
      cat(paste("tau = ", tau[i], "\n", sep = ""))
      printCoefmat(x$tTable[[i]], signif.stars = TRUE, P.values = TRUE)
      cat("\n")
    }
  }

  invisible(x)
}
