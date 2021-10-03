ranef.lqmm <- function(object, ...) {
  tau <- object$tau
  nq <- length(tau)
  group <- object$group
  M <- object$ngroups
  BLPu <- vector("list", M)
  q <- object$dim_theta[2]
  mmr.l <- split(object$mmr, group)
  sigma <- VarCorr(object)
  cov_name <- object$cov_name

  if (cov_name %in% c("pdIdent", "pdDiag")) {
    if (nq == 1) {
      sigma <- diag(x = sigma, nrow = length(sigma))
    } else {
      for (j in 1:nq) sigma[[j]] <- diag(x = sigma[[j]], nrow = length(sigma[[j]]))
    }
  }

  if (nq == 1) {
    psi <- varAL(object$scale, tau)
    INV <- lapply(mmr.l, function(x, a, b, q) {
      x <- matrix(x, ncol = q)
      n <- nrow(x)
      y <- x %*% a %*% t(x) + diag(b, n)
      solve(y)
    }, a = sigma, b = psi, q = q)
    GZ <- lapply(mmr.l, function(x, a, q) {
      x <- matrix(x, ncol = q)
      a %*% t(x)
    }, a = sigma, q = q)
    RES <- split(object$y - object$mmf %*% matrix(object$theta_x) - meanAL(0, object$scale, tau), group)
    for (i in 1:M) {
      BLPu[[i]] <- GZ[[i]] %*% INV[[i]] %*% matrix(RES[[i]])
    }
    ans <- data.frame(matrix(unlist(BLPu), ncol = q, byrow = TRUE))
    rownames(ans) <- unique(group)
    colnames(ans) <- object$mm
  } else {
    ans <- vector("list", nq)
    for (j in 1:nq) {
      tmp <- object[[j]]
      psi <- varAL(tmp$scale, tau[j])
      INV <- lapply(mmr.l, function(x, a, b, q) {
        x <- matrix(x, ncol = q)
        n <- nrow(x)
        y <- x %*% a %*% t(x) + diag(b, n)
        solve(y)
      },
      a = sigma[[j]], b = psi, q = q
      )
      GZ <- lapply(mmr.l, function(x, a, q) {
        x <- matrix(x, ncol = q)
        a %*% t(x)
      }, a = sigma[[j]], q = q)
      RES <- split(object$y - object$mmf %*% matrix(tmp$theta_x) - meanAL(0, tmp$scale, tau[j]), group)
      for (i in 1:M) {
        BLPu[[i]] <- GZ[[i]] %*% INV[[i]] %*% matrix(RES[[i]])
      }
      ans[[j]] <- data.frame(matrix(unlist(BLPu), ncol = q, byrow = TRUE))
      rownames(ans[[j]]) <- unique(group)
      colnames(ans[[j]]) <- object$mm
    }
    names(ans) <- format(tau, digits = 4)
  }
  return(ans)
}
