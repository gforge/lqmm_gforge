
covHandling <- function(theta, n, cov_name, quad_type) {
  if (cov_name %in% c("pdIdent", "pdDiag")) {
    if (quad_type == "robust") {
      sigma <- theta
      if (any(sigma < 0)) {
        warning("Not positive-definite variance-covariance of random effects.")
        sigma[sigma < 0] <- .Machine$double.eps
      }
      sigma <- varAL(sigma, 0.5)
    } else {
      sigma <- theta
      if (any(sigma < 0)) {
        warning("Not positive-definite variance-covariance of random effects.")
        sigma[sigma < 0] <- .Machine$double.eps
      }
      sigma <- sigma^2
    }
  }


  if (cov_name == "pdCompSymm") {
    if (quad_type == "robust") {
      stop("Not implemented yet: Gauss-Laguerre quadrature requires uncorrelated random effects.")
    } else {
      sigma <- as.matrix(invTfun(x = theta, n = n, type = cov_name))
      sigma <- sigma %*% sigma
      if (!is.positive.definite(sigma)) {
        warning("Not positive-definite variance-covariance of random effects.")
        sigma <- make.positive.definite(sigma)
      }
    }
  }

  if (cov_name == "pdSymm") {
    if (quad_type == "robust") {
      stop("Not implemented yet: Gauss-Laguerre quadrature requires uncorrelated random effects.")
    } else {
      sigma <- as.matrix(invTfun(x = theta, n = n, type = cov_name))
      sigma <- sigma %*% sigma
      if (!is.positive.definite(sigma)) {
        warning("Not positive-definite variance-covariance of random effects.")
        sigma <- make.positive.definite(sigma)
      }
    }
  }

  return(sigma)
}
