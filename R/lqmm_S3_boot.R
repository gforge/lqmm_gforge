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
