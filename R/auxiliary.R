###            Fit linear quantile models and linear quantile mixed models
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

.onAttach <- function(lib, pkg) {
  if (interactive() || getOption("verbose")) {
    packageStartupMessage(sprintf(
      "Package %s (%s) loaded. Type citation(\"%s\") on how to cite this package\n", pkg,
      packageDescription(pkg)$Version, pkg
    ))
  }
}


#' Retrieves an expression from parent environments
#'
#' As we update models objects may end up in different environments
#' and we need to be able to locate them regardless of where they exist.
#' This function traverses through all parent environments in search for
#' the object of interest.
#'
#' @param expr The expression to search for
#'
#' @return The evaluated object
retrieveExpression <- function(expr) {
  for (level in 1:sys.nframe()) {
    out <- tryCatch(eval(expr, envir = parent.frame(n = level)), error = \(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  stop("Failed to retrieve value")
}

##

allVarsRec <- function(object) {

  # source `nlme' version 3.1-100 (Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the R Development Core Team)

  if (is.list(object)) {
    unlist(lapply(object, allVarsRec))
  } else {
    all.vars(object)
  }
}

asOneFormula <- function(..., omit = c(".", "pi")) {
  # source `nlme' version 3.1-100 (Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the R Development Core Team)

  names <- unique(allVarsRec((list(...))))
  names <- names[is.na(match(names, omit))]
  if (length(names)) {
    eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  } else {
    NULL
  }
}

##

bandwidth.rq <- function(p, n, hs = TRUE, alpha = 0.05) {

  # source `quantreg' version 5.11 (Roger Koenker)

  x0 <- qnorm(p)
  f0 <- dnorm(x0)
  if (hs) {
    n^(-1 / 3) * qnorm(1 - alpha / 2)^(2 / 3) * ((1.5 * f0^2) / (2 *
      x0^2 + 1))^(1 / 3)
  } else {
    n^-0.2 * ((4.5 * f0^4) / (2 * x0^2 + 1)^2)^0.2
  }
}

##########################################################################################
# Asymmetric Laplace distribution



#' The Asymmetric Laplace Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the asymmetric Laplace distribution.
#'
#'
#' The asymmetric Laplace distribution with parameters (mu, sigma, tau) has
#' density \deqn{f(x) = \tau(1-\tau)/\sigma e^{-1/(2\sigma) (\theta max(x,0) +
#' (1 - \theta) max(-x,0))}}
#'
#' @aliases pal qal ral dal
#' @param x vector of quantiles (\code{dal}, \code{pal}) or probabilities
#' (\code{qal}).
#' @param n number of observations.
#' @param mu location parameter.
#' @param sigma positive scale parameter.
#' @param tau skewness parameter (0,1).
#' @param log logical; if \code{TRUE}, probabilities are log--transformed.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}, \code{\link{lqm}}
#' @keywords asymmetric Laplace distribution
dal <- function(x, mu = 0, sigma = 1, tau = 0.5, log = FALSE) {
  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (tau == 0) tau <- eps
  if (tau == 1) tau <- 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")

  ind <- ifelse(x < mu, 1, 0)

  val <- tau * (1 - tau) / sigma * exp(-(x - mu) / sigma * (tau - ind))

  if (log) log(val) else val
}

pal <- function(x, mu = 0, sigma = 1, tau = 0.5) {
  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (tau == 0) tau <- eps
  if (tau == 1) tau <- 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")

  ifelse(x < mu, tau * exp((1 - tau) / sigma * (x - mu)),
    1 - (1 - tau) * exp(-tau / sigma * (x - mu))
  )
}

ral <- function(n, mu = 0, sigma = 1, tau = 0.5) {
  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (tau == 0) tau <- eps
  if (tau == 1) tau <- 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")

  u <- runif(n)

  x1 <- mu + sigma / (1 - tau) * log(u / tau)

  x2 <- mu - sigma / tau * log((1 - u) / (1 - tau))

  ifelse(x1 < mu, x1, x2)
}

qal <- function(x, mu = 0, sigma = 1, tau = 0.5) {
  if (any(x > 1) | any(x < 0)) stop("x must be in [0,1]")

  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (tau == 0) tau <- eps
  if (tau == 1) tau <- 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")

  ifelse(x < tau, mu + (sigma / (1 - tau)) * log(x / tau),
    mu - (sigma / tau) * log((1 - x) / (1 - tau))
  )
}



#' Functions for Asymmetric Laplace Distribution Parameters
#'
#' Accessory functions.
#'
#' \code{meanAL} computes the mean of an asymmetric Laplace with parameters
#' \code{mu}, \code{sigma} and \code{tau}.
#'
#' \code{varAL} computes the variance of an asymmetric Laplace with parameters
#' \code{sigma} and \code{tau}.
#'
#' \code{invvarAL} computes the scale parameter of an asymmetric Laplace with
#' parameter \code{tau} and variance \code{x}.
#'
#' @aliases varAL invvarAL meanAL
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param tau skewness parameter.
#' @param x numeric value.
#' @author Marco Geraci
#' @seealso \code{\link{dal}}, \code{\link{mleAL}}
#' @references Yu K and Zhang J (2005). A three-parameter asymmetric Laplace
#' distribution and its extension. Communications in Statistics-Theory and
#' Methods 34, 1867--1879.
#' @keywords asymmetric Laplace distribution maximum likelihood estimation
meanAL <- function(mu, sigma, tau) {
  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (tau == 0) tau <- eps
  if (tau == 1) tau <- 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")

  mu + sigma * (1 - 2 * tau) / (tau * (1 - tau))
}

varAL <- function(sigma, tau) {
  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (any(tau == 0)) tau[tau == 0] <- eps
  if (any(tau == 1)) tau[tau == 1] <- 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")

  sigma^2 * (1 - 2 * tau + 2 * tau^2) / ((1 - tau)^2 * tau^2)
}

invvarAL <- function(x, tau) {
  eps <- .Machine$double.eps^(2 / 3)
  if (any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
  if (tau == 0) tau <- eps
  if (tau == 1) tau <- 1 - eps

  sqrt(x * (tau * (1 - tau))^2 / (1 - 2 * tau + 2 * tau^2))
}



#' Maximum Likelihood Estimation of Asymmetric Laplace Distribution
#'
#' This function estimates the parameters of an asymmetric Laplace distribution
#' for a sample.
#'
#'
#' @param x a numeric vector.
#' @return
#'
#' an object of class \code{list} containing the following components:
#'
#' \item{m}{location parameter} \item{sigma}{scale parameter}
#' \item{tau}{skewness parameter} \item{r}{number of iterations}
#' @author Marco Geraci
#' @seealso \code{\link{dal}}, \code{\link{meanAL}}
#' @references Yu K and Zhang J (2005). A three-parameter asymmetric Laplace
#' distribution and its extension. Communications in Statistics-Theory and
#' Methods 34, 1867--1879.
#' @keywords asymmetric Laplace distribution maximum likelihood estimation
mleAL <- function(x) {
  tau <- 0.5
  m <- as.numeric(quantile(x, tau))
  sigma <- 1

  r <- 0
  while (r < 1000) {
    m.last <- as.numeric(quantile(x, tau))
    res <- x - m.last
    sigma.last <- mean(res * (tau - ifelse(res < 0, 1, 0)))
    a <- mean(res * ifelse(res <= 0, 1, 0))
    tau.last <- (a + sqrt(a^2 - (mean(x) - m.last) * a)) / (mean(x) - m.last)

    dm <- abs(m - m.last)
    ds <- abs(sigma - sigma.last)
    dp <- abs(tau - tau.last)

    if (all(c(dm, ds, dp) < 0.0001)) {
      break
    } else {
      m <- m.last
      tau <- tau.last
      sigma <- sigma.last
    }
    r <- r + 1
  }
  list(m = m, sigma = sigma, tau = tau, r = r)
}
