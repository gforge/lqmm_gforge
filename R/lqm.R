###            Fit a linear quantile model (continuous and count reponses)
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


##########################################################################################
# lqm functions (independent data)
# (negative) Laplace log-likelihood, gradient and Hessian

lqm <- function(formula, data, subset, na.action, weights = NULL, tau = 0.5, contrasts = NULL, control = list(), fit = TRUE) {
  if (any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
  nq <- length(tau)

  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) {
    stop("'weights' must be a numeric vector")
  }
  if (is.null(w)) w <- rep(1, length(y))
  x <- model.matrix(mt, mf, contrasts)
  dim_theta <- ncol(x)

  # Control

  if (is.null(names(control))) {
    control <- lqmControl()
  } else {
    control_default <- lqmControl()
    control_names <- intersect(names(control), names(control_default))
    control_default[control_names] <- control[control_names]
    control <- control_default
  }
  if (is.null(control$loop_step)) control$loop_step <- sd(as.numeric(y))
  if (control$beta > 1 || control$beta < 0) stop("Beta must be a decreasing factor in (0,1)")
  if (control$gamma < 1) stop("Beta must be a nondecreasing factor >= 1")


  # Starting values

  theta_0 <- lm.wfit(x = as.matrix(x), y = y, w = w)$coefficients

  if (!fit) {
    return(list(theta = theta_0, x = as.matrix(x), y = y, weights = w, tau = tau, control = control))
  }

  if (nq == 1) {
    fit <- lqm.fit.gs(theta = theta_0, x = as.matrix(x), y = y, weights = w, tau = tau, control = control)
  } else {
    fit <- vector("list", nq)
    names(fit) <- format(tau, digits = 4)
    for (i in 1:nq) fit[[i]] <- lqm.fit.gs(theta = theta_0, x = as.matrix(x), y = y, weights = w, tau = tau[i], control = control)
  }

  term.labels <- colnames(x)

  if (nq > 1) {
    fit$theta <- matrix(NA, dim_theta, nq)
    fit$scale <- rep(NA, nq)

    for (i in 1:nq) {
      fit$theta[, i] <- fit[[i]]$theta
      fit$scale[i] <- fit[[i]]$scale
    }
    rownames(fit$theta) <- term.labels
    colnames(fit$theta) <- format(tau, digits = 4)
  }

  class(fit) <- "lqm"

  fit$call <- Call
  fit$na.action <- attr(mf, "na.action")
  fit$contrasts <- attr(x, "contrasts")
  fit$term.labels <- term.labels
  fit$terms <- mt
  fit$nobs <- length(y)
  fit$dim_theta <- fit$edf <- dim_theta
  fit$rdf <- fit$nobs - fit$edf
  fit$tau <- tau
  fit$x <- as.matrix(x)
  fit$y <- y
  fit$weights <- w
  fit$levels <- .getXlevels(mt, mf)
  fit$InitialPar <- list(theta = theta_0)
  fit$control <- control

  fit
}
