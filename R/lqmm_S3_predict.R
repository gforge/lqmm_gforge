#' Predictions from an \code{lqmm} Object
#'
#' The predictions at level 0 correspond to predictions based only on the fixed
#' effects estimates. The predictions at level 1 are obtained by adding the
#' best linear predictions of the random effects to the predictions at level 0.
#' See details for interpretation. The function \code{predint} will produce
#' 1-alpha confidence intervals based on bootstrap centiles.
#'
#'
#' As discussed by Geraci and Bottai (2014), integrating over the random
#' effects will give "weighted averages" of the cluster-specific quantile
#' effects. These may be interpreted strictly as population regression
#' quantiles only for the median (\code{tau=0.5}). Therefore, predictions at
#' the population level (\code{code=0}) should be interpreted analogously.
#'
#' @aliases predint predint.lqmm predict.lqmm
#' @param object an \code{lqmm} object.
#' @param level an optional integer vector giving the level of grouping to be
#' used in obtaining the predictions.
#' @param alpha 1-\code{alpha} is the confidence level.
#' @param R number of bootstrap replications.
#' @param seed optional random number generator seed.
#' @param \dots not used.
#' @return a vector or a matrix of predictions for \code{predict.lqmm}. A data
#' frame or a list of data frames for \code{predint.lqmm} containing
#' predictions, lower and upper bounds of prediction intervals, and standard
#' errors.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}, \code{\link{ranef.lqmm}},
#' \code{\link{coef.lqmm}}
#' @references Geraci M and Bottai M (2014). Linear quantile mixed models.
#' Statistics and Computing, 24(3), 461--479.
#' @keywords prediction
#' @examples
#'
#' ## Orthodont data
#' data(Orthodont)
#'
#' # Random intercept model
#' fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
#' 	tau = c(0.1,0.5,0.9), data = Orthodont)
#'
#' # Predict (y - Xb)
#' predict(fitOi.lqmm, level = 0)
#'
#' # Predict (y - Xb - Zu)
#' predict(fitOi.lqmm, level = 1)
#'
#' # 95% confidence intervals
#' predint(fitOi.lqmm, level = 0, alpha = 0.05)
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
