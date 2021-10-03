Tfun <- function(n, type = "pdSymm") {
  val <- 0

  if (type == "pdIdent") {
    val <- matrix(diag(n), ncol = 1)
  }

  if (type == "pdDiag") {
    val <- sapply(1:n, function(x, n) {
      z <- matrix(0, n, n)
      z[x, x] <- 1
      z
    }, n = n)
  }

  if (type == "pdSymm") {
    val <- apply(diag(n * (n + 1) / 2), 2, invTfun, n = n, type = type)
  }

  if (type == "pdCompSymm") {
    A <- matrix(1, n, n)
    diag(A) <- rep(0, n)
    val <- if (n > 1) cbind(as.vector(diag(n)), as.vector(A)) else 1
  }

  return(val)
}

invTfun <- function(x, n, type = "pdSymm") {
  val <- NULL

  if (type == "pdCompSymm") {
    val <- matrix(x[2], n, n)
    diag(val) <- rep(x[1], n)
  }

  if (type == "pdSymm") {
    dim(x) <- NULL
    val <- matrix(0, n, n)
    val[lower.tri(val, diag = TRUE)] <- x
    hold <- val
    hold[upper.tri(hold, diag = TRUE)] <- 0
    val <- val + t(hold)
  }

  return(val)
}
