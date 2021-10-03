#' Gaussian Quadrature
#'
#' This function calculates nodes and weights for Gaussian quadrature. See
#' \code{help("gauss.quad")} from package \code{statmod}.
#'
#'
#' @param n number of nodes and weigh
#' @param kind kind of Gaussian quadrature, one of "legendre", "chebyshev1", "chebyshev2",
#'  "hermite", "jacobi" or "laguerre"
#' @param alpha parameter for Jacobi or Laguerre quadrature, must be greater than -1
#' @param beta parameter for Jacobi quadrature, must be greater than -1
#'
#' @author Original version by Gordon Smyth
#' @source Gordon Smyth with contributions from Yifang Hu, Peter Dunn and
#' Belinda Phipson. (2011). statmod: Statistical Modeling. R package version
#' 1.4.11. \url{https://CRAN.R-project.org/package=statmod}
#' @keywords gaussian quadrature
gauss.quad <- function(n, kind = "legendre", alpha = 0, beta = 0) {

  # source `statmod' version 1.4.11 (Gordon Smyth)

  n <- as.integer(n)
  if (n < 0) {
    stop("need non-negative number of nodes")
  }
  if (n == 0) {
    return(list(nodes = numeric(0), weights = numeric(0)))
  }
  kind <- match.arg(kind, c(
    "legendre", "chebyshev1", "chebyshev2",
    "hermite", "jacobi", "laguerre"
  ))
  i <- 1:n
  i1 <- i[-n]
  switch(kind,
    legendre = {
      muzero <- 2
      a <- rep(0, n)
      b <- i1 / sqrt(4 * i1^2 - 1)
    },
    chebyshev1 = {
      muzero <- pi
      a <- rep(0, n)
      b <- rep(0.5, n - 1)
      b[1] <- sqrt(0.5)
    },
    chebyshev2 = {
      muzero <- pi / 2
      a <- rep(0, n)
      b <- rep(0.5, n - 1)
    },
    hermite = {
      muzero <- sqrt(pi)
      a <- rep(0, n)
      b <- sqrt(i1 / 2)
    },
    jacobi = {
      ab <- alpha + beta
      muzero <- 2^(ab + 1) * gamma(alpha + 1) * gamma(beta +
        1) / gamma(ab + 2)
      a <- i
      a[1] <- (beta - alpha) / (ab + 2)
      i2 <- 2:n
      abi <- ab + 2 * i2
      a[i2] <- (beta^2 - alpha^2) / (abi - 2) / abi
      b <- i1
      b[1] <- sqrt(4 * (alpha + 1) * (beta + 1) / (ab + 2)^2 / (ab +
        3))
      i2 <- i1[-1]
      abi <- ab + 2 * i2
      b[i2] <- sqrt(4 * i2 * (i2 + alpha) * (i2 + beta) * (i2 +
        ab) / (abi^2 - 1) / abi^2)
    },
    laguerre = {
      a <- 2 * i - 1 + alpha
      b <- sqrt(i1 * (i1 + alpha))
      muzero <- gamma(alpha + 1)
    }
  )
  A <- rep(0, n * n)
  A[(n + 1) * (i - 1) + 1] <- a
  A[(n + 1) * (i1 - 1) + 2] <- b
  A[(n + 1) * i1] <- b
  dim(A) <- c(n, n)
  vd <- eigen(A, symmetric = TRUE)
  w <- rev(as.vector(vd$vectors[1, ]))
  w <- muzero * w^2
  x <- rev(vd$values)
  list(nodes = x, weights = w)
}



#' Gaussian Quadrature
#'
#' This function calculates nodes and weights for Gaussian quadrature in terms
#' of probability distributions. See \code{help("gauss.quad.prob")} from
#' package \code{statmod}.
#'
#' @param dist distribution that Gaussian quadrature is based on, one of "uniform", "normal",
#'  "beta" or "gamma"
#' @param l lower limit of uniform distribution
#' @param u upper limit of uniform distribution
#' @param mu mean of normal distribution
#' @param sigma standard deviation of normal distribution
#' @param alpha positive shape parameter for gamma distribution or first shape parameter for beta
#'  distribution
#' @param beta positive scale parameter for gamma distribution or second shape parameter for
#'  beta distribution
#'
#' @author Original version by Gordon Smyth
#' @source Gordon Smyth with contributions from Yifang Hu, Peter Dunn and
#' Belinda Phipson. (2011). statmod: Statistical Modeling. R package version
#' 1.4.11. \url{https://CRAN.R-project.org/package=statmod}
#' @keywords gaussian quadrature
#' @inheritParams gauss.quad
gauss.quad.prob <- function(n, dist = "uniform", l = 0, u = 1, mu = 0, sigma = 1,
                              alpha = 1, beta = 1) {

  # source `statmod' version 1.4.11 (Gordon Smyth)

  n <- as.integer(n)
  if (n < 0) {
    stop("need non-negative number of nodes")
  }
  if (n == 0) {
    return(list(nodes = numeric(0), weights = numeric(0)))
  }
  dist <- match.arg(dist, c(
    "uniform", "beta1", "beta2", "normal",
    "beta", "gamma"
  ))
  if (n == 1) {
    switch(dist,
      uniform = {
        x <- (l + u) / 2
      },
      beta1 = ,
      beta2 = ,
      beta = {
        x <- alpha / (alpha + beta)
      },
      normal = {
        x <- mu
      },
      gamma = {
        x <- alpha * beta
      }
    )
    return(list(nodes = x, weights = 1))
  }
  if (dist == "beta" && alpha == 0.5 && beta == 0.5) {
    dist <- "beta1"
  }
  if (dist == "beta" && alpha == 1.5 && beta == 1.5) {
    dist <- "beta2"
  }
  i <- 1:n
  i1 <- 1:(n - 1)
  switch(dist,
    uniform = {
      a <- rep(0, n)
      b <- i1 / sqrt(4 * i1^2 - 1)
    },
    beta1 = {
      a <- rep(0, n)
      b <- rep(0.5, n - 1)
      b[1] <- sqrt(0.5)
    },
    beta2 = {
      a <- rep(0, n)
      b <- rep(0.5, n - 1)
    },
    normal = {
      a <- rep(0, n)
      b <- sqrt(i1 / 2)
    },
    beta = {
      ab <- alpha + beta
      a <- i
      a[1] <- (alpha - beta) / ab
      i2 <- 2:n
      abi <- ab - 2 + 2 * i2
      a[i2] <- ((alpha - 1)^2 - (beta - 1)^2) / (abi - 2) / abi
      b <- i1
      b[1] <- sqrt(4 * alpha * beta / ab^2 / (ab + 1))
      i2 <- i1[-1]
      abi <- ab - 2 + 2 * i2
      b[i2] <- sqrt(4 * i2 * (i2 + alpha - 1) * (i2 + beta -
        1) * (i2 + ab - 2) / (abi^2 - 1) / abi^2)
    },
    gamma = {
      a <- 2 * i + alpha - 2
      b <- sqrt(i1 * (i1 + alpha - 1))
    }
  )
  A <- rep(0, n * n)
  A[(n + 1) * (i - 1) + 1] <- a
  A[(n + 1) * (i1 - 1) + 2] <- b
  A[(n + 1) * i1] <- b
  dim(A) <- c(n, n)
  vd <- eigen(A, symmetric = TRUE)
  w <- rev(as.vector(vd$vectors[1, ]))^2
  x <- rev(vd$values)
  switch(dist,
    uniform = x <- l + (u - l) * (x + 1) / 2,
    beta1 = ,
    beta2 = ,
    beta = x <- (x + 1) / 2,
    normal = x <- mu + sqrt(2) *
      sigma * x,
    gamma = x <- beta * x
  )
  list(nodes = x, weights = w)
}
