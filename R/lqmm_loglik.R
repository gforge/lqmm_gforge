loglik.t <- function(theta, sigma, x, y, z, weights, Tq, V, W, tau, p, q, m, M, N, Kq, minn, maxn) {
  if (length(theta) != (p + m)) stop("Check length theta")

  ans <- .C("C_ll_h", theta = as.double(theta), as.double(x), as.double(y), as.double(z), as.double(weights), as.double(Tq), as.double(V), as.double(W), as.double(sigma), as.single(tau), as.integer(p), as.integer(q), as.integer(m), as.integer(M), as.integer(N), as.integer(Kq), as.integer(minn - 1), as.integer(maxn), loglik = double(1)) # , PACKAGE = "lqmm")

  ans$loglik
}

loglik.s <- function(sigma, theta, x, y, z, weights, Tq, V, W, tau, p, q, m, M, N, Kq, minn, maxn) {
  if (length(theta) != (p + m)) stop("Check length theta")

  ans <- .C("C_ll_h", theta = as.double(theta), as.double(x), as.double(y), as.double(z), as.double(weights), as.double(Tq), as.double(V), as.double(W), as.double(sigma), as.single(tau), as.integer(p), as.integer(q), as.integer(m), as.integer(M), as.integer(N), as.integer(Kq), as.integer(minn - 1), as.integer(maxn), loglik = double(1)) # , PACKAGE = "lqmm")

  ans$loglik
}
