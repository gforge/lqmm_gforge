logliki <- function(theta, x, y, tau, smooth = FALSE, omicron = 0.001) {
  n <- length(y)
  res <- y - x %*% matrix(theta)

  if (smooth) {
    s <- ifelse(res <= (tau - 1) * omicron, -1, ifelse(res >= tau * omicron, 1, 0))
    w <- as.numeric(1 - s^2)
    W <- diag(w, n, n)
    gs <- s * ((2 * tau - 1) * s + 1) / 2
    cs <- sum(0.25 * (1 - 2 * tau) * omicron * s - 0.25 * (1 - 2 * tau + 2 * tau^2) * omicron * s^2)
    res <- matrix(res)
    ans <- as.numeric(1 / (2 * omicron) * t(res) %*% W %*% res + t(gs) %*% res + cs)
    grad <- -t(x) %*% matrix(1 / omicron * W %*% res + gs)
    hess <- 1 / omicron * t(x) %*% W %*% x
  } else {
    ind <- tau - as.numeric(res < 0)
    ans <- as.numeric(sum(res * ind))
    grad <- -t(x) %*% ind
    hess <- matrix(0, ncol(x), ncol(x))
  }

  if (is.na(ans)) {
    ans <- Inf
    grad <- matrix(Inf, ncol(x))
  }

  attr(ans, "grad") <- matrix(grad)

  return(ans)
}
