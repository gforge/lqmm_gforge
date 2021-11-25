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
#' @param object an \code{lqmm} object.
#' @param level an optional integer vector giving the level of grouping to be
#' used in obtaining the predictions.
#' @param alpha 1-\code{alpha} is the confidence level.
#' @param R number of bootstrap replications.
#' @param seed optional random number generator seed.
#' @param newdata Defaults to the original data used for generating the model.
#'  Currently only implemented for `level = 0` predictions.
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
#' @rdname lqmm_predictions
#' @examples
#'
#' ## Orthodont data
#' data(Orthodont)
#'
#' # Random intercept model
#' fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject, tau = c(0.1,0.5,0.9), data = Orthodont)
#'
#' # Predict (y - Xb)
#' predict(fitOi.lqmm, level = 0)
#'
#' # Predict (y - Xb - Zu)
#' predict(fitOi.lqmm, level = 1)
#'
#' # 95% confidence intervals
#' predint(fitOi.lqmm, level = 0, alpha = 0.05)
predict.lqmm <- function(object,
                         newdata,
                         level = 0,
                         na.action = na.pass,
                         ...)
{
  tau <- object$tau
  nq <- length(tau)
  q <- object$dim_theta[2]
  if (!level %in% c(0, 1)) stop("level must be either 0 (population-averaged) or 1 (conditional)")

  if (!missing(newdata)) {
    ## check newdata
    if (!inherits(newdata, "data.frame")) stop("'newdata' must be a data frame")
    if (!all(all.vars(object$mtf) %in% names(newdata))) stop("newdata must have all terms in 'fixed' formula from main call")
    if (!all(all.vars(object$mtr) %in% names(newdata))) stop("newdata must have all terms in 'random' formula from main call")

    ## ordering data by groups
    if (level == 1) {
      groupFormula <- asOneSidedFormula(attr(object$group, "name"))
      grp <- model.frame(groupFormula, newdata)
      origOrder <- row.names(newdata)
      ord <- order(unlist(grp, use.names = FALSE))
      grp <- grp[ord,,drop = TRUE]
      newdata <- newdata[ord, ,drop = FALSE]
      revOrder <- match(origOrder, row.names(newdata)) # putting in orig. order
    } else {
      grp <- rep(1, times = nrow(newdata))
      revOrder <- 1:nrow(newdata)
    }

    ## create data frames
    mtf <- object$mtf
    mtr <- object$mtr
    mtf <- delete.response(mtf)
    mf <- model.frame(formula(mtf), newdata, na.action = na.action, drop.unused.levels = TRUE, xlev = object$xlevels[['fixed']])
    mr <- model.frame(formula(mtr), newdata, na.action = na.action, drop.unused.levels = TRUE, xlev = object$xlevels[['random']])

    if (!is.null(cl <- attr(mtf, "dataClasses")))
      .checkMFClasses(cl, mf)
    if (!is.null(cl <- attr(mtr, "dataClasses")))
      .checkMFClasses(cl, mr)
    object$mmf <- model.matrix(formula(mtf), mf)
    object$mmr <- model.matrix(formula(mtr), mr)

    object$group <- grp
    object$revOrder <- revOrder
  }
  group <- object$group
  M <- object$ngroups


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
        u <- as.numeric(RE.l[[match(names(mmr.l)[i], names(RE.l))]])
        RND <- rbind(RND, matrix(as.numeric(mmr.l[[i]]), ncol = q) %*% matrix(u, nrow = q))
      } # for i
    } else {
      RND <- matrix(NA, length(object$y), nq)
      for (j in 1:nq) {
        RE.l <- split(RE[[j]], unique(group))
        tmp <- NULL
        for (i in 1:M) {
          u <- as.numeric(RE.l[[match(names(mmr.l)[i], names(RE.l))]])
          tmp <- rbind(tmp, matrix(as.numeric(mmr.l[[i]]), ncol = q) %*% matrix(u, nrow = q))
        } # for i
        RND[,j] <- tmp
      } # for j
    } # else nq
  } # if level

  if (level == 0) {
    colnames(FXD) <- format(tau, digits = 4)
    ans <- FXD[object$revOrder,]
  }

  if (level == 1) {
    ans <- FXD + RND
    colnames(ans) <- format(tau, digits = 4)
    ans <- ans[object$revOrder,]
  }

  return(ans)
}
