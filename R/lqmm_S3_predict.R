predict.lqmm <- function(object, level = 0, ...) {
  tau <- object$tau
  nq <- length(tau)
  group <- object$group
  M <- object$ngroups
  q <- object$dim_theta[2]

  if (nq == 1) {
    FXD <- object$mmf %*% matrix(object$theta_x)
  } else {
    FXD <- object$mmf %*% object$theta_x
  }

  if (level == 1) {
    RE <- ranef(object)
    mmr.l <- split(object$mmr, group)
    if (nq == 1) {
      RE.l <- split(RE, unique(group))
      RND <- NULL
      for (i in 1:M) {
        RND <- rbind(RND, matrix(as.numeric(mmr.l[[i]]), ncol = q) %*% matrix(as.numeric(RE.l[[i]]), nrow = q))
      }
    } else {
      RND <- matrix(NA, length(object$y), nq)
      for (j in 1:nq) {
        RE.l <- split(RE[[j]], unique(group))
        tmp <- NULL
        for (i in 1:M) {
          tmp <- rbind(tmp, matrix(as.numeric(mmr.l[[i]]), ncol = q) %*% matrix(as.numeric(RE.l[[i]]), nrow = q))
        }
        RND[, j] <- tmp
      }
    }
  }

  if (level == 0) {
    colnames(FXD) <- format(tau, digits = 4)
    ans <- FXD[object$revOrder, ]
  }
  if (level == 1) {
    ans <- FXD + RND
    colnames(ans) <- format(tau, digits = 4)
    ans <- ans[object$revOrder, ]
  }

  return(ans)
}
