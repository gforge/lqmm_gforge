#' Extract Random Effects
#'
#' This function computes random effects for a linear quantile mixed model.
#'
#' The prediction of the random effects is done via estimated best linear
#' prediction (Geraci and Bottai, 2014). The generic function \code{ranef} is
#' imported from the \code{nlme} package (Pinheiro et al, 2014).
#'
#' @aliases ranef ranef.lqmm
#' @param object an object of \code{\link{class}} \code{lqmm}.
#' @param \dots not used.
#' @return a data frame or a list of data frames of predicted random effects.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}, \code{\link{coef.lqmm}}
#' @references Geraci M and Bottai M (2014). Linear quantile mixed models.
#' Statistics and Computing, 24(3), 461--479. doi: 10.1007/s11222-013-9381-9.
#'
#' Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear
#' and Nonlinear Mixed Effects Models. R package version 3.1-117,
#' \url{https://CRAN.R-project.org/package=nlme}.
#' @keywords random effects coefficients
#' @export
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
