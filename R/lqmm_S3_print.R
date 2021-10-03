print.lqmm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  tau <- x$tau
  nq <- length(tau)

  if (nq == 1) {
    theta_x <- x$theta_x
    names(theta_x) <- x$nn
    sigma <- VarCorr(x)
    psi <- varAL(x$scale, tau)

    cat("Call: ")
    dput(x$call)
    cat("\n")
    cat(paste("Quantile", tau, "\n"))
    cat("\n")

    cat("Fixed effects:\n")
    print.default(format(theta_x, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Covariance matrix of the random effects:\n")
    print.default(format(sigma, digits = digits), quote = FALSE)

    cat("\n")
    cat(paste("Residual scale parameter: ", format(x$scale, digits = digits),
      " (standard deviation ", format(sqrt(psi), digits = digits), ")", "\n",
      sep = ""
    ))
    cat(paste("Log-likelihood:", format(x$logLik, digits = digits), "\n"))
    cat(paste("\nNumber of observations:", length(x$y), "\n"))
    cat(paste("Number of groups:", x$ngroups, "\n"))
  } else {
    theta_x <- x$theta_x
    colnames(theta_x) <- paste("tau = ", format(tau, digits = digits), sep = "")
    rownames(theta_x) <- x$nn
    Scale <- sapply(x[1:nq], function(z) z$scale)
    psi <- varAL(sigma = Scale, tau = tau)
    sigma <- VarCorr(x)
    ll <- sapply(x[1:nq], function(z) z$logLik)

    cat("Call: ")
    dput(x$call)
    cat("\n")
    cat("Fixed effects:\n")
    print.default(format(theta_x, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Covariance matrix of the random effects:\n")
    for (i in 1:nq) {
      cat(paste("tau = "), tau[i], "\n", sep = "")
      print.default(format(sigma[[i]], digits = digits), quote = FALSE)
    }

    cat("\n")
    cat("Residual scale parameter: ")
    cat(paste(format(Scale, digits = digits), " (tau = ", tau, ") ", sep = ""))
    cat("\n")
    cat("Log-likelihood: ")
    cat(paste(format(ll, digits = digits), " (tau = ", tau, ") ", sep = ""))
    cat("\n")
    cat(paste("\nNumber of observations:", length(x$y), "\n"))
    cat(paste("Number of groups:", x$ngroups, "\n"))
  }

  invisible(x)
}
