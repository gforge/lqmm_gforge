lqm.fit.gs <- function(theta, x, y, weights, tau, control) {
  n <- length(y)
  p <- ncol(x)
  wx <- x * weights
  wy <- y * weights

  if (is.null(p)) stop("x must be a matrix")
  if (missing(theta)) theta <- lm.fit(as.matrix(wx), wy)$coefficients
  if (is.null(control$loop_step)) control$loop_step <- sd(as.numeric(wy))
  if (length(tau) > 1) {
    tau <- tau[1]
    warning("Length of tau is greater than 1. Only first value taken")
  }

  if (control$method == "gs1") {
    fit <- .C("C_gradientSi", theta = as.double(theta), as.double(wx), as.double(wy), as.single(tau), as.integer(n), as.integer(p), as.double(control$loop_step), as.double(control$beta), as.double(control$gamma), as.integer(control$reset_step), as.double(control$loop_tol_ll), as.double(control$loop_tol_theta), as.integer(control$check_theta), as.integer(control$loop_max_iter), as.integer(control$verbose), CONVERGE = integer(1), grad = double(p), optimum = double(1)) # , PACKAGE = "lqmm")
  } else if (control$method == "gs2") {
    fit <- gradientSi(theta, wx, wy, tau, control)
  }

  fit$residuals <- y - x %*% matrix(fit$theta)
  fit$scale <- weighted.mean(fit$residuals * (tau - (fit$residuals < 0)), weights)
  fit$logLik <- n * log(tau * (1 - tau) / fit$scale) - 1 / fit$scale * fit$optimum
  OPTIMIZATION <- list(loop = fit$CONVERGE)

  errorHandling(OPTIMIZATION$loop, "low", control$loop_max_iter, control$loop_tol_ll, "lqm")

  list(theta = fit$theta, scale = fit$scale, gradient = matrix(fit$grad), logLik = fit$logLik, opt = OPTIMIZATION)
}
