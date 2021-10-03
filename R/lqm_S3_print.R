print.lqm <- function(x, digits = max(6, getOption("digits")), ...) {
  tau <- x$tau
  nq <- length(tau)

  if (nq == 1) {
    theta <- x$theta
    names(theta) <- x$term.labels
    psi <- varAL(x$scale, tau)

    cat("Call: ")
    dput(x$call)
    cat("\n")
    cat(paste("Quantile", tau, "\n"))
    cat("Fixed effects:\n")
    print.default(format(theta, digits = digits), print.gap = 2, quote = FALSE)

    cat("\nDegrees of freedom:", x$nobs, "total;", x$rdf, "residual\n")
    cat(paste("Residual scale parameter: ", format(x$scale, digits = digits),
      " (standard deviation ", format(sqrt(psi), digits = digits), ")", "\n",
      sep = ""
    ))
    cat(paste("Log-likelihood (Laplace):", format(x$logLik, digits = digits), "\n"))
  } else {
    theta <- x$theta
    rownames(theta) <- x$term.labels
    colnames(theta) <- paste("tau = ", format(tau, digits = digits), sep = "")

    cat("Call: ")
    dput(x$call)
    cat("\n")
    cat("Fixed effects:\n")
    print.default(format(theta, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nDegrees of freedom:", x$nobs, "total;", x$rdf, "residual\n")
  }

  invisible(x)
}
