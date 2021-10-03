gradientSi <- function(theta, x, y, tau, control) {
  step <- control$loop_step
  maxiter <- control$loop_max_iter
  omicron <- control$omicron

  theta_0 <- theta
  ll_0 <- logliki(theta_0, x, y, tau, control$smooth, omicron)
  eps <- .Machine$double.eps

  for (i in 1:maxiter) {
    if (control$verbose) cat(paste0("  (", i, ") logLik = ", round(ll_0, 12), "\n"))
    # line search
    theta_1 <- theta_0 - attributes(ll_0)$grad * step
    ll_1 <- logliki(theta_1, x, y, tau, control$smooth, omicron)
    if (ll_1 > ll_0) {
      if (control$verbose) cat("  Decreasing step...\n")
      step <- step * control$beta
    } else {
      rule <- if (control$check_theta) "2" else "1"
      check <- switch_check(ll_0, ll_1, control$loop_tol_ll, theta_0, theta_1, control$loop_tol_theta, rule = rule)
      if (check) break
      theta_0 <- theta_1
      ll_0 <- ll_1
      step <- if (control$reset_step) control$loop_step else step * control$gamma
    }
  }

  list(theta = as.numeric(theta_1), grad = attributes(ll_1)$grad, optimum = as.numeric(ll_1), CONVERGE = if (i == maxiter) -1 else i)
}
