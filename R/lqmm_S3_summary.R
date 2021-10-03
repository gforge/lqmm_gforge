#' Summary for an \code{lqmm} Object
#'
#' Summary method for class \code{lqmm}.
#'
#' \code{print.summary.lqmm} formats the coefficients, standard errors, etc.
#' and additionally gives `significance stars'.
#'
#' @param object an object of \code{\link{class}} \code{lqmm}.
#' @param method specifies the method used to compute standard errors.
#' Currently, only the bootstrap method ("boot") is available.
#' @param alpha significance level.
#' @param covariance logical flag. If \code{TRUE} the covariance matrix of the
#' bootstrap estimates is provided.
#' @param \dots see \code{\link{boot.lqmm}} for additional arguments.
#' @return an object of class \code{summary.lqmm}. The function
#' \code{summary.lqmm} computes and returns a list of summary statistics of the
#' fitted linear quantile mixed model given in \code{object}, using the
#' components (list elements) from its argument, plus
#'
#' \item{Cov}{the covariance matrix obtained from the bootstrapped estimates
#' (if \code{covariance = TRUE}).} \item{tTable}{a matrix with estimates,
#' standard errors, etc.} \item{B}{the matrix of all bootstrapped parameters.}
#' @author Marco Geraci
#' @seealso \code{\link{print.summary.lqmm}} \code{\link{lqmm}}
#' @keywords bootstrap standard errors
#' @examples
#'
#' data(Orthodont)
#' fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
#' 	tau = c(0.1,0.5,0.9), data = Orthodont)
#' summary(fitOi.lqmm)
summary.lqmm <- function(object, method = "boot", alpha = 0.05, covariance = FALSE, ...) {
  tau <- object$tau
  nq <- length(tau)
  object$logLik <- logLik(object)
  object$aic <- AIC(object)
  est <- extractAll(object)
  nn <- object$nn
  rdf <- object$rdf
  ddf <- object$dim_theta[1] - 1
  npars <- object$dim_theta[1] + object$dim_theta_z + 1

  coln <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")

  if (ddf > 0) {
    newF <- update.formula(as.formula(object$call[["fixed"]]), as.formula(. ~ 1))
    FITNULL <- update(object, fixed = as.formula(newF), evaluate = TRUE)
    LR <- -2 * (logLik(FITNULL) - object$logLik)
    if (any(LR < 0)) {
      LR[LR < 0] <- 0
      warning("Negative LR test (set to zero).")
    }
    LRp <- 1 - pchisq(LR, df = ddf)
    names(LR) <- tau
    attr(LR, "pvalue") <- LRp
    attr(LR, "df") <- ddf
    object$LRtest <- LR
    object$null_model <- FITNULL
  }

  if (method == "boot") {
    B <- boot.lqmm(object, ...)
    R <- attr(B, "R")
    if (nq == 1) {
      Cov <- cov(as.matrix(B))
      stds <- sqrt(diag(Cov))
      tP <- 2 * pt(-abs(est / stds), R - 1)
      lower <- est + qt(alpha / 2, R - 1) * stds
      upper <- est + qt(1 - alpha / 2, R - 1) * stds
      ans <- cbind(est, stds, lower, upper, tP)
      colnames(ans) <- coln
      ans <- ans[c(nn), , drop = FALSE]
    } else {
      Cov <- apply(B, 3, function(x) cov(as.matrix(x)))
      stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
      tP <- 2 * pt(-abs(est / stds), R - 1)
      lower <- est + qt(alpha / 2, R - 1) * stds
      upper <- est + qt(1 - alpha / 2, R - 1) * stds
      ans <- vector("list", nq)
      names(ans) <- tau
      Cov.array <- array(NA, dim = c(object$df, object$df, nq))
      for (i in 1:nq) {
        ans[[i]] <- cbind(est[, i], stds[, i], lower[, i], upper[, i], tP[, i])
        rownames(ans[[i]]) <- rownames(est)
        colnames(ans[[i]]) <- coln
        ans[[i]] <- ans[[i]][c(nn), , drop = FALSE]
        Cov.array[, , i] <- matrix(Cov[, i], object$df, object$df)
      }
      Cov <- Cov.array
      dimnames(Cov) <- list(rownames(attr(B, "estimated")), rownames(attr(B, "estimated")), format(tau, digits = 4))
    }
  }
  if (method == "nid") {



  }
  if (covariance) object$Cov <- Cov
  object$tTable <- ans
  object$B <- B
  class(object) <- "summary.lqmm"
  return(object)
}

#' Print an \code{lqmm} Summary Object
#'
#' Print summary of an \code{lqmm} object.
#'
#'
#' @param x a \code{summary.lqmm} object.
#' @param digits a non-null value for digits specifies the minimum number of
#' significant digits to be printed in values.
#' @param \dots not used.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}, \code{\link{summary.lqmm}}
#' @keywords print summary
print.summary.lqmm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  tau <- x$tau
  nq <- length(tau)

  cat("Call: ")
  dput(x$call)
  cat("\n")

  if (nq == 1) {
    cat(paste("Quantile", tau, "\n"))
    cat("\n")

    cat("Fixed effects:\n")
    printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
  } else {
    for (i in 1:nq) {
      cat(paste("tau = ", tau[i], "\n", sep = ""))
      cat("\n")
      cat("Fixed effects:\n")

      printCoefmat(x$tTable[[i]], signif.stars = TRUE, P.values = TRUE)
      cat("\n")
    }
  }

  cat("AIC:\n")
  print.default(
    paste(format(x$aic, digits = digits), " (df = ", attr(x$logLik, "df"), ")", sep = ""),
    quote = FALSE
  )
}
