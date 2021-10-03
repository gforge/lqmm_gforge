score.al <- function(theta, x, y, tau, smooth, omicron = 0.001) {
  res <- as.numeric(y - x %*% matrix(theta))
  n <- length(y)

  if (smooth) {
    s <- ifelse(res <= (tau - 1) * omicron, -1, ifelse(res >= tau * omicron, 1, 0))
    w <- as.numeric(1 - s^2)
    W <- diag(w, n, n)
    gs <- s * ((2 * tau - 1) * s + 1) / 2
    ans <- -t(x) %*% matrix(1 / omicron * W %*% res + gs)
  } else {
    ans <- -t(x) %*% matrix(tau - as.numeric(res < 0), n, 1)
  }

  if (any(is.na(ans))) {
    ans <- matrix(Inf, ncol(x))
  }


  return(matrix(ans))
}
