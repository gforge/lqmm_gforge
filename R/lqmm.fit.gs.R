#' Linear Quantile Mixed Models Fitting by Gradient Search
#'
#' This function controls the arguments to be passed to routines written in C
#' for LQMM estimation. The optimization algorithm is based on the gradient of
#' the Laplace log--likelihood (Bottai, Orsini and Geraci, 2014; Geraci and
#' Bottai, 2014).
#'
#' In \code{\link{lqmm}}, see argument \code{fit} for generating a list of
#' arguments to be called by this function; see argument \code{covariance} for
#' alternative variance--covariance matrices.
#'
#' NOTE: the data should be ordered by \code{group} when passed to
#' \code{lqmm.fit.gs} (such ordering is performed by \code{\link{lqmm}}).
#'
#' @param theta_0 starting values for the linear predictor.
#' @param x the model matrix for fixed effects (see details).
#' @param y the model response (see details).
#' @param z the model matrix for random effects (see details).
#' @param weights the weights used in the fitting process (see details).
#' @param cov_name variance--covariance matrix of the random effects. Default
#' is \code{pdIdent}. See details.
#' @param V nodes of the quadrature.
#' @param W weights of the quadrature.
#' @param sigma_0 starting value for the scale parameter.
#' @param tau the quantile(s) to be estimated.
#' @param group the grouping factor (see details).
#' @param control list of control parameters used for optimization (see
#' \code{\link{lqmmControl}}).
#' @return An object of class "list" containing the following components:
#'
#' \item{theta}{a vector of coefficients, including the "raw"
#' variance--covariance parameters (see \code{\link{VarCorr.lqmm}}).}
#' \item{scale}{the scale parameter.} \item{gradient}{the gradient.}
#' \item{logLik}{the log--likelihood.} \item{opt}{number of iterations when the
#' estimation algorithm stopped for lower (theta) and upper (scale) loop.}.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}
#' @references Bottai M, Orsini N, Geraci M. (2014). A gradient search
#' maximization algorithm for the asymmetric Laplace likelihood, Journal of
#' Statistical Computation and Simulation (in press).
#'
#' Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and
#' Computing, 24(3), 461--479.
#' @keywords fitting
#' @examples
#'
#' set.seed(123)
#'
#' M <- 50
#' n <- 10
#' test <- data.frame(x = runif(n*M,0,1), group = rep(1:M,each=n))
#' test$y <- 10*test$x + rep(rnorm(M, 0, 2), each = n) + rchisq(n*M, 3)
#' lqmm.ls <- lqmm(fixed = y ~ x, random = ~ 1, group = group,
#' 	data = test, fit = FALSE)
#'
#' do.call("lqmm.fit.gs", lqmm.ls)
#'
#'
lqmm.fit.gs <- function(theta_0, x, y, z, weights, cov_name, V, W, sigma_0, tau, group, control) {
  if (length(tau) > 1) {
    tau <- tau[1]
    warning("Length of tau is greater than 1. Only first value taken")
  }

  x <- as.matrix(x)
  z <- as.matrix(z)
  V <- as.matrix(V)
  W <- as.matrix(W)

  p <- ncol(x)
  q <- ncol(z)
  m <- theta.z.dim(cov_name, q)
  # cov_type <- cov.sel(cov_name)
  Tq <- Tfun(n = q, type = cov_name)

  N <- nrow(x)
  if (length(y) != N) stop("Check dim of x and y")
  Kq <- nrow(V)

  ns <- as.integer(table(group))
  M <- length(ns)
  minn <- c(1, cumsum(ns[-M]) + 1)
  maxn <- cumsum(ns)

  if (length(weights) != M) stop("Length of \"weights\" does not match number of groups")

  if (is.null(control$LP_step)) control$LP_step <- sd(as.numeric(y))
  UP_max_iter <- control$UP_max_iter
  if (UP_max_iter == 0) stop("Increase number of maximum iterations", " (", UP_max_iter, ")", sep = "")
  r <- 0

  while (r < UP_max_iter) {
    if (control$verbose) cat(paste("Upper loop = ", r + 1, "\n", sep = ""))

    ans <- .C("C_gradientSh",
      theta = as.double(theta_0), as.double(x), as.double(y), as.double(z), as.double(weights), as.double(Tq), as.double(V), as.double(W),
      as.double(sigma_0), as.single(tau), as.integer(p), as.integer(q), as.integer(m), as.integer(M), as.integer(N), as.integer(Kq), as.integer(minn - 1),
      as.integer(maxn), as.double(control$LP_step), as.double(control$beta), as.double(control$gamma), as.integer(control$reset_step),
      as.double(control$LP_tol_ll), as.double(control$LP_tol_theta), as.integer(control$check_theta), as.integer(control$LP_max_iter),
      as.integer(control$verbose), low_loop = integer(1), double(1), grad = double(p + m), opt_val = double(1)
    ) # , PACKAGE = "lqmmGforge")

    theta_1 <- ans$theta
    grad <- ans$grad

    opt_s <- optimize(f = loglik.s, interval = c(.Machine$double.eps, 10 * sigma_0), theta = theta_1, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn)

    sigma_1 <- opt_s$minimum

    delta <- abs(sigma_1 - sigma_0)

    if (delta < control$UP_tol) {
      break
    } else {
      r <- r + 1
      theta_0 <- theta_1
      sigma_0 <- sigma_1
    }
  }

  if (r < UP_max_iter) upp_loop <- r + 1
  if (r == UP_max_iter & UP_max_iter > 0) upp_loop <- -1
  if (r == UP_max_iter & UP_max_iter == 0) upp_loop <- -2
  OPTIMIZATION <- list(low_loop = ans$low_loop, upp_loop = upp_loop)

  errorHandling(OPTIMIZATION$low_loop, "low", control$LP_max_iter, control$LP_tol_ll, "lqmm")
  errorHandling(OPTIMIZATION$upp_loop, "upp", control$UP_max_iter, control$UP_tol, "lqmm")

  list(theta = theta_1, scale = sigma_1, gradient = grad, logLik = -ans$opt_val, opt = OPTIMIZATION)
}
