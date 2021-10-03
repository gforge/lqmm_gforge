lqmm.fit.df <- function(theta_0, x, y, z, weights, cov_name, V, W, sigma_0, tau, group, control) {
  if (length(tau) > 1) {
    tau <- tau[1]
    warning("Length of tau is greater than 1. Only first value taken")
  }

  x <- as.matrix(x)
  z <- as.matrix(z)
  V <- as.matrix(V)
  W <- as.matrix(W)

  p <- ncol(x)
  q <- ncol(z)
  m <- theta.z.dim(cov_name, q)
  # cov_type <- cov.sel(cov_name)
  Tq <- Tfun(n = q, type = cov_name)

  N <- nrow(x)
  if (length(y) != N) stop("Check dim of x and y")
  Kq <- nrow(V)

  ns <- as.integer(table(group))
  M <- length(ns)
  minn <- c(1, cumsum(ns[-M]) + 1)
  maxn <- cumsum(ns)

  if (length(weights) != M) stop("Length of \"weights\" does not match number of groups")

  UP_max_iter <- control$UP_max_iter
  if (UP_max_iter == 0) stop("Increase number of maximum iterations", " (", UP_max_iter, ")", sep = "")
  r <- 0

  while (r < UP_max_iter) {
    if (control$verbose) cat(paste("Upper loop = ", r + 1, "\n", sep = ""))

    ans <- optim(par = theta_0, fn = loglik.t, sigma = sigma_0, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn, control = list(maxit = control$LP_max_iter, abstol = control$LP_tol_ll, trace = control$verbose))

    if (control$verbose) cat(paste("(", r + 1, ")", " logLik = ", round(-ans$value, 3), "\n", sep = ""))

    theta_1 <- ans$par

    opt_s <- optimize(f = loglik.s, interval = c(.Machine$double.eps, 1e3 * sigma_0), theta = theta_1, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn)

    sigma_1 <- opt_s$minimum

    delta <- abs(sigma_1 - sigma_0)

    if (delta < control$UP_tol) {
      break
    } else {
      r <- r + 1
      theta_0 <- theta_1
      sigma_0 <- sigma_1
    }
  }

  low_loop <- ans$convergence

  if (r < UP_max_iter) upp_loop <- r + 1
  if (r == UP_max_iter & UP_max_iter > 0) upp_loop <- -1
  if (r == UP_max_iter & UP_max_iter == 0) upp_loop <- -2
  low_loop <- if (low_loop == 1) -1 else as.numeric(ans$counts[1])
  OPTIMIZATION <- list(low_loop = low_loop, upp_loop = upp_loop)

  errorHandling(OPTIMIZATION$low_loop, "low", control$LP_max_iter, control$LP_tol_ll, "lqmm")
  errorHandling(OPTIMIZATION$upp_loop, "upp", control$UP_max_iter, control$UP_tol, "lqmm")

  list(theta = theta_1, scale = sigma_1, logLik = -ans$value, opt = OPTIMIZATION)
}
