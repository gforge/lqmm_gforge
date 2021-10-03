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



#' Fitting Linear Quantile Models
#' 
#' \code{lqm} is used to fit linear quantile models based on the asymmetric
#' Laplace distribution.
#' 
#' The function computes an estimate on the tau-th quantile function of the
#' response, conditional on the covariates, as specified by the formula
#' argument. The quantile predictor is assumed to be linear. The function
#' maximizes the (log)likelihood of a Laplace regression which is equivalent to
#' the minimization of the weighted sum of absolute residuals (Koenker and
#' Bassett, 1978). The optimization algorithm is based on the gradient of the
#' Laplace log--likelihood (Bottai, Orsini and Geraci, 2013).
#' 
#' @param formula an object of class \code{\link{formula}} for fixed effects: a
#' symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in
#' the model. If not found in data, the variables are taken from
#' \code{environment(formula)}, typically the environment from which \code{lqm}
#' is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is set by the \code{na.action} setting of
#' \code{\link{options}}.
#' @param weights An optional vector of weights to be used in the fitting
#' process.
#' @param tau the quantile(s) to be estimated. This must be a number between 0
#' and 1, otherwise the execution is stopped. If more than one quantile is
#' specified, rounding off to the 4th decimal must give non--duplicated values
#' of \code{tau}, otherwise the execution is stopped.
#' @param contrasts an optional list. See the contrasts.arg of
#' \code{\link{model.matrix.default}}.
#' @param control list of control parameters of the fitting process. See
#' \code{\link{lqmControl}}.
#' @param fit logical flag. If \code{FALSE} the function returns a list of
#' arguments to be passed to \code{lqm.fit.gs}.
#' @return \code{lqm} returns an object of \code{\link{class}} \code{lqm}.
#' 
#' The function \code{summary} is used to obtain and print a summary of the
#' results.
#' 
#' An object of class \code{lqm} is a list containing the following components:
#' 
#' \item{theta}{a vector of coefficients. \code{theta} is a named matrix of
#' coefficients when \code{tau} is a vector of values.} \item{scale}{the scale
#' parameter.} \item{gradient}{the gradient.} \item{logLik}{the
#' log--likelihood.} \item{opt}{details on optimization (see
#' \code{\link{lqm.fit.gs}}).} \item{call}{the matched call.}
#' \item{term.labels}{names for theta.} \item{terms}{the terms object used.}
#' \item{nobs}{the number of observations.} \item{edf,dim_theta}{the length of
#' theta.} \item{rdf}{the number of residual degrees of freedom.}
#' \item{tau}{the estimated quantile(s).} \item{x}{the model matrix.}
#' \item{y}{the model response.} \item{weights}{the weights used in the fitting
#' process (a vector of 1's if \code{weights} = NULL).}
#' \item{InitialPar}{starting values for theta.} \item{control}{list of control
#' parameters used for optimization (see \code{\link{lqmControl}}).}
#' @note Updates/FAQ/news are published here
#' \url{http://marcogeraci.wordpress.com/}. New versions are usually published
#' here \url{https://r-forge.r-project.org/R/?group_id=1396} before going on
#' CRAN.
#' @author Marco Geraci
#' @seealso \code{\link{summary.lqm}, \link{coef.lqm}, \link{predict.lqm},
#' \link{residuals.lqm}}
#' @references Bottai M, Orsini N, Geraci M (2015). A Gradient Search
#' Maximization Algorithm for the Asymmetric Laplace Likelihood, Journal of
#' Statistical Computation and Simulation, 85(10), 1919-1925.
#' 
#' Chen C (2007). A finite smoothing algorithm for quantile regression. Journal
#' of Computational and Graphical Statistics, 16(1), 136-164.
#' 
#' Koenker R and Bassett G (1978). Regression Quantiles. Econometrica 46(1),
#' 33--50.
#' @keywords quantile regression
#' @examples
#' 
#' 
#' set.seed(123)
#' n <- 500
#' p <- 1:3/4
#' test <- data.frame(x = runif(n,0,1))
#' test$y <- 30 + test$x + rnorm(n)
#' fit.lqm <- lqm(y ~ x, data = test, tau = p,
#' 	control = list(verbose = FALSE, loop_tol_ll = 1e-9), fit = TRUE)
#' fit.lqm
#' 
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
