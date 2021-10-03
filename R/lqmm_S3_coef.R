#' Extract LQMM Coefficients
#' 
#' \code{coef} extracts model coefficients from \code{lqmm} objects.
#' 
#' 
#' @param object a fitted object of \code{\link{class}} "lqmm".
#' @param \dots not used.
#' @return a vector for single quantiles or a matrix for multiple quantiles.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}} \code{\link{summary.lqmm}}
#' @keywords coefficients
coef.lqmm <- function(object, ...) {
  tau <- object$tau
  nq <- length(tau)
  ans <- object$theta_x

  if (nq == 1) {
    names(ans) <- object$nn
  }

  return(ans)
}



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
#' 
#' 
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

predint.lqmm <- function(object, level = 0, alpha = 0.05, R = 50, seed = round(runif(1, 1, 10000))) {
  tau <- object$tau
  nq <- length(object$tau)
  p <- object$dim_theta[1]
  m <- object$dim_theta_z

  B <- boot(object, R = R, seed = seed)
  tmp <- object

  if (nq == 1) {
    yhat <- matrix(NA, object$nobs, R)
    for (i in 1:R) {
      tmp$theta <- B[i, 1:(p + m)]
      tmp$theta_x <- B[i, 1:p]
      tmp$theta_z <- B[i, (p + 1):(p + m)]
      tmp$scale <- B[i, (p + m + 1)]
      yhat[, i] <- predict(tmp, level = level)
    }
    LB <- apply(yhat, 1, quantile, probs = alpha / 2)
    UB <- apply(yhat, 1, quantile, probs = 1 - alpha / 2)
    ans <- data.frame(yhat = predict(object, level = level), lower = LB, upper = UB, SE = apply(yhat, 1, sd))
  } else {
    ans <- list()
    for (j in 1:nq) {
      tmp$tau <- tau[j]
      yhat <- matrix(NA, object$nobs, R)
      for (i in 1:R) {
        tmp$theta <- B[i, 1:(p + m), j]
        tmp$theta_x <- B[i, 1:p, j]
        tmp$theta_z <- B[i, (p + 1):(p + m), j]
        tmp$scale <- B[i, (p + m + 1), j]
        yhat[, i] <- predict(tmp, level = level)
      }
      LB <- apply(yhat, 1, quantile, probs = alpha / 2)
      UB <- apply(yhat, 1, quantile, probs = 1 - alpha / 2)
      ans[[j]] <- data.frame(yhat = predict(tmp, level = level), lower = LB, upper = UB, SE = apply(yhat, 1, sd))
    }
    names(ans) <- format(tau, digits = 4)
  }

  return(ans)
}



#' Residuals from an \code{lqmm} Object
#' 
#' The residuals at level 0 correspond to population residuals (based only on
#' the fixed effects estimates). The residuals at level 1 are obtained by
#' adding the best linear predictions of the random effects to the predictions
#' at level 0 and the subtracting these from the model response.
#' 
#' 
#' @param object an \code{lqmm} object.
#' @param level an optional integer vector giving the level of grouping to be
#' used in obtaining the predictions. Level zero corresponds to the population
#' residuals.
#' @param \dots not used.
#' @return a matrix of residuals.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}, \code{\link{predict.lqmm}},
#' \code{\link{coef.lqmm}}, \code{\link{ranef.lqmm}},
#' @references Geraci M and Bottai M (2014). Linear quantile mixed models.
#' Statistics and Computing, 24(3), 461--479. doi: 10.1007/s11222-013-9381-9.
#' @keywords residuals
residuals.lqmm <- function(object, level = 0, ...) {
  object$y[object$revOrder] - predict(object, level = level)
}



#' Extract Log-Likelihood
#' 
#' \code{logLik.lqmm} extracts the log-likelihood of a fitted LQMM.
#' 
#' 
#' @param object an object of \code{\link{class}} "lqmm".
#' @param \dots not used.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}} \code{\link{AIC}}
#' @keywords models
logLik.lqmm <- function(object, ...) {
  tdf <- object$edf + 1
  tau <- object$tau
  nq <- length(tau)

  if (nq == 1) {
    ans <- object$logLik
  } else {
    ans <- NULL
    for (i in 1:nq) ans <- c(ans, object[[i]]$logLik)
    names(ans) <- as.character(format(tau, digits = 4))
  }

  attr(ans, "nobs") <- object$nobs
  attr(ans, "df") <- tdf
  attr(ans, "class") <- "logLik"

  return(ans)
}



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
#' 
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

boot.lqmm <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE) {
  if (startQR) warning("Standard errors may be underestimated when 'startQR = TRUE'")

  set.seed(seed)
  tau <- object$tau
  nq <- length(tau)

  ngroups <- object$ngroups
  group_all <- object$group
  group_unique <- unique(group_all)
  obsS <- replicate(R, sample(group_unique, replace = TRUE))
  dim_theta_z <- object$dim_theta_z

  npars <- object$dim_theta[1] + dim_theta_z + 1
  dimn <- c(object$nn, paste("reStruct", 1:dim_theta_z, sep = ""), "scale")

  control <- object$control
  control$verbose <- FALSE

  FIT_ARGS <- list(
    x = as.matrix(object$mmf), y = object$y, z = as.matrix(object$mmr), cov_name = object$cov_name,
    V = object$QUAD$nodes, W = object$QUAD$weights, group = group_all, control = control
  )

  if (nq == 1) {
    FIT_ARGS$theta_0 <- object$theta
    FIT_ARGS$sigma_0 <- object$scale
    FIT_ARGS$tau <- object$tau

    bootmat <- matrix(NA, R, npars)
    colnames(bootmat) <- dimn

    for (i in 1:R) {
      group_freq <- table(obsS[, i])
      sel_unique <- group_unique %in% names(group_freq)
      w <- rep(0, ngroups)
      w[sel_unique] <- group_freq
      FIT_ARGS$weights <- w

      if (!startQR) {
        lmfit <- lm.wfit(x = as.matrix(object$mmf), y = object$y, w = rep(w, table(group_all)))
        theta_x <- lmfit$coefficients
        theta_z <- if (object$type == "normal") {
          rep(1, dim_theta_z)
        } else {
          rep(invvarAL(1, 0.5), dim_theta_z)
        }
        FIT_ARGS$theta_0 <- c(theta_x, theta_z)
        FIT_ARGS$sigma_0 <- invvarAL(mean(lmfit$residuals^2), 0.5)
      }

      if (control$method == "gs") fit <- try(do.call(lqmm.fit.gs, FIT_ARGS))
      if (control$method == "df") fit <- try(do.call(lqmm.fit.df, FIT_ARGS))
      if (!inherits(fit, "try-error")) bootmat[i, ] <- c(fit$theta, fit$scale)
    }
  } else {
    bootmat <- array(NA, dim = c(R, npars, nq), dimnames = list(NULL, dimn, paste("tau = ", format(tau, digits = 4), sep = "")))

    for (i in 1:R) {
      group_freq <- table(obsS[, i])
      sel_unique <- group_unique %in% names(group_freq)
      w <- rep(0, ngroups)
      w[sel_unique] <- group_freq
      FIT_ARGS$weights <- w
      for (j in 1:nq) {
        if (startQR) {
          FIT_ARGS$theta_0 <- object[[j]]$theta
          FIT_ARGS$sigma_0 <- object[[j]]$scale
        } else {
          lmfit <- lm.wfit(x = as.matrix(object$mmf), y = object$y, w = rep(w, table(group_all)))
          theta_x <- lmfit$coefficients
          theta_z <- if (object$type == "normal") {
            rep(1, dim_theta_z)
          } else {
            rep(invvarAL(1, 0.5), dim_theta_z)
          }
          FIT_ARGS$theta_0 <- c(theta_x, theta_z)
          FIT_ARGS$sigma_0 <- invvarAL(mean(lmfit$residuals^2), 0.5)
        }

        FIT_ARGS$tau <- object$tau[j]

        if (control$method == "gs") fit <- try(do.call(lqmm.fit.gs, FIT_ARGS))
        if (control$method == "df") fit <- try(do.call(lqmm.fit.df, FIT_ARGS))
        if (!inherits(fit, "try-error")) bootmat[i, , j] <- c(fit$theta, fit$scale)
      }
    }
  }

  class(bootmat) <- "boot.lqmm"
  attr(bootmat, "tau") <- tau
  attr(bootmat, "estimated") <- extractAll(object)
  attr(bootmat, "R") <- R
  attr(bootmat, "seed") <- seed
  attr(bootmat, "nn") <- object$nn
  attr(bootmat, "npars") <- npars
  attr(bootmat, "indices") <- obsS
  attr(bootmat, "dim_theta") <- object$dim_theta
  attr(bootmat, "dim_theta_z") <- object$dim_theta_z

  return(bootmat)
}

extractBoot.boot.lqmm <- function(object, which = "fixed") {
  tau <- attr(object, "tau")
  nq <- length(tau)
  nn <- attr(object, "nn")
  dim_theta <- attr(object, "dim_theta")
  dim_theta_z <- attr(object, "dim_theta_z")

  ind.f <- 1:dim_theta[1]
  ind.r <- (dim_theta[1] + 1):(dim_theta[1] + dim_theta_z)
  ind.s <- dim_theta[1] + dim_theta_z + 1

  if (which == "fixed") {
    if (nq == 1) {
      ans <- as.matrix(object[, c(ind.f, ind.s)])
    } else {
      ans <- object[, c(ind.f, ind.s), ]
    }
  }

  if (which == "random") {
    if (nq == 1) {
      ans <- as.matrix(object[, ind.r])
    } else {
      ans <- object[, ind.r, ]
    }
  }

  return(ans)
}

summary.boot.lqmm <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), ...) {
  est <- attr(object, "estimated")
  R <- attr(object, "R")
  tau <- attr(object, "tau")
  nq <- length(tau)
  nn <- attr(object, "nn")
  npars <- attr(object, "npars")
  coln <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")

  if (nq == 1) {
    Cov <- cov(as.matrix(object))
    stds <- sqrt(diag(Cov))
    tP <- 2 * pt(-abs(est / stds), R - 1)
    lower <- est + qt(alpha / 2, R - 1) * stds
    upper <- est + qt(1 - alpha / 2, R - 1) * stds
    ans <- cbind(est, stds, lower, upper, tP)
    colnames(ans) <- coln
    ans <- ans[c(nn, "scale"), ]
    cat(paste("Quantile", tau, "\n"))
    printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
  } else {
    Cov <- apply(object, 3, function(x) cov(as.matrix(x)))
    stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
    tP <- 2 * pt(-abs(est / stds), R - 1)
    lower <- est + qt(alpha / 2, R - 1) * stds
    upper <- est + qt(1 - alpha / 2, R - 1) * stds
    for (i in 1:nq) {
      ans <- cbind(est[, i], stds[, i], lower[, i], upper[, i], tP[, i])
      rownames(ans) <- rownames(est)
      colnames(ans) <- coln
      ans <- ans[c(nn, "scale"), ]
      cat(paste("tau = ", tau[i], "\n", sep = ""))
      printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
      cat("\n")
    }
  }
}

extractAll <- function(object) {
  nq <- length(object$tau)

  dim_theta_z <- object$dim_theta_z
  dimn <- c(object$nn, paste("reStruct", 1:dim_theta_z, sep = ""), "scale")

  if (nq == 1) {
    ans <- c(object$theta, object$scale)
    names(ans) <- dimn
  } else {
    val <- NULL
    for (i in 1:nq) val <- c(val, object[[i]]$scale)
    ans <- rbind(object$theta_x, object$theta_z, val)
    rownames(ans) <- dimn
  }
  return(ans)
}
