VarCorr.lqmm <- function(x, sigma = NULL, ...) {
  tau <- x$tau
  nq <- length(tau)

  theta_z <- x$theta_z
  dim_theta <- x$dim_theta
  q <- dim_theta[2]
  cov_name <- x$cov_name
  type <- x$type
  mm <- x$mm

  if (nq == 1) {
    sigma <- covHandling(theta = theta_z, n = q, cov_name = cov_name, quad_type = type)
    if (cov_name == "pdIdent") {
      sigma <- rep(sigma, q)
      names(sigma) <- mm
    }
    if (cov_name == "pdDiag") {
      names(sigma) <- mm
    }
    if (cov_name == "pdCompSymm") {
      rownames(sigma) <- colnames(sigma) <- mm
    }
    if (cov_name == "pdSymm") {
      rownames(sigma) <- colnames(sigma) <- mm
    }
  } else {
    sigma <- vector("list", nq)
    names(sigma) <- format(tau, digits = 4)
    for (i in 1:nq) {
      sigma[[i]] <- covHandling(theta = theta_z[, i], n = q, cov_name = cov_name, quad_type = type)
      if (cov_name == "pdIdent") {
        sigma[[i]] <- rep(sigma[[i]], q)
        names(sigma[[i]]) <- mm
      }
      if (cov_name == "pdDiag") {
        names(sigma[[i]]) <- mm
      }
      if (cov_name == "pdCompSymm") {
        rownames(sigma[[i]]) <- colnames(sigma[[i]]) <- mm
      }
      if (cov_name == "pdSymm") {
        rownames(sigma[[i]]) <- colnames(sigma[[i]]) <- mm
      }
    }
  }

  return(sigma)
}
