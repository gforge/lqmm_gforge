#' Linear Quantile Mixed Models Fitting by Derivative-Free Optimization
#' 
#' This function controls the arguments to be passed to \code{\link{optim}} and
#' \code{\link{optimize}} for LQMM estimation.
#' 
#' In \code{\link{lqmm}}, see argument \code{fit} for generating a list of
#' arguments to be called by this function; see argument \code{covariance} for
#' alternative variance--covariance matrices.
#' 
#' NOTE: the data should be ordered by \code{group} when passed to
#' \code{lqmm.fit.df} (such ordering is performed by \code{\link{lqmm}}).
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
#' \item{scale}{the scale parameter.} \item{logLik}{the log--likelihood.}
#' \item{opt}{number of iterations when the estimation algorithm stopped for
#' lower (theta) and upper (scale) loop.}.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}
#' @keywords fitting
#' @examples
#' 
#' set.seed(123)
#' 
#' M <- 50
#' n <- 10
#' test <- data.frame(x = runif(n*M,0,1), group = rep(1:M,each=n))
#' test$y <- 10*test$x + rep(rnorm(M, 0, 2), each = n) + rchisq(n*M, 3)
#' lqmm.ls <- lqmm(fixed = y ~ x, random = ~ 1, group = group, data = test,
#' 	fit = FALSE)
#' 
#' do.call("lqmm.fit.df", lqmm.ls)
#' 
#' 
lqmm.fit.df <- function(theta_0, x, y, z, weights, cov_name, V, W, sigma_0, tau, group, control) {
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

  UP_max_iter <- control$UP_max_iter
  if (UP_max_iter == 0) stop("Increase number of maximum iterations", " (", UP_max_iter, ")", sep = "")
  r <- 0

  while (r < UP_max_iter) {
    if (control$verbose) cat(paste("Upper loop = ", r + 1, "\n", sep = ""))

    ans <- optim(par = theta_0, fn = loglik.t, sigma = sigma_0, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn, control = list(maxit = control$LP_max_iter, abstol = control$LP_tol_ll, trace = control$verbose))

    if (control$verbose) cat(paste("(", r + 1, ")", " logLik = ", round(-ans$value, 3), "\n", sep = ""))

    theta_1 <- ans$par

    opt_s <- optimize(f = loglik.s, interval = c(.Machine$double.eps, 1e3 * sigma_0), theta = theta_1, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn)

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

  low_loop <- ans$convergence

  if (r < UP_max_iter) upp_loop <- r + 1
  if (r == UP_max_iter & UP_max_iter > 0) upp_loop <- -1
  if (r == UP_max_iter & UP_max_iter == 0) upp_loop <- -2
  low_loop <- if (low_loop == 1) -1 else as.numeric(ans$counts[1])
  OPTIMIZATION <- list(low_loop = low_loop, upp_loop = upp_loop)

  errorHandling(OPTIMIZATION$low_loop, "low", control$LP_max_iter, control$LP_tol_ll, "lqmm")
  errorHandling(OPTIMIZATION$upp_loop, "upp", control$UP_max_iter, control$UP_tol, "lqmm")

  list(theta = theta_1, scale = sigma_1, logLik = -ans$value, opt = OPTIMIZATION)
}
