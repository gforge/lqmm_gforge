#' Control parameters for lqmm estimation
#' 
#' A list of parameters for controlling the fitting process.
#' 
#' \code{LP} (lower loop) refers to the estimation of regression coefficients
#' and variance-covariance parameters. \code{UP} (upper loop) refers to the
#' estimation of the scale parameter.
#' 
#' @param method character vector that specifies the estimation method: "gs"
#' for gradient search (default) and "df" for Nelder-Mead.
#' @param LP_tol_ll tolerance expressed as absolute change of the
#' log-likelihood.
#' @param LP_tol_theta tolerance expressed as absolute change of \code{theta}
#' @param check_theta logical flag. If TRUE the algorithm performs an
#' additional check on the change in the estimates.
#' @param LP_step step size (default standard deviation of response).
#' @param beta decreasing step factor for line search (0,1).
#' @param gamma nondecreasing step factor for line search (>= 1).
#' @param reset_step logical flag. If \code{TRUE} the step size is reset to the
#' initial value at each iteration.
#' @param LP_max_iter maximum number of iterations
#' @param UP_tol tolerance expressed as absolute change of the \code{scale}
#' parameter.
#' @param UP_max_iter maximum number of iterations.
#' @param startQR logical flag. If \code{FALSE} (default) the least squares
#' estimate of the fixed effects is used as starting value of \code{theta_x}
#' and \code{scale}. If \code{TRUE} the \code{\link{lqm}} estimate is used.
#' @param verbose logical flag.
#' @return a list of control parameters.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}}
#' @keywords fitting control
lqmmControl <- function(method = "gs", LP_tol_ll = 1e-5, LP_tol_theta = 1e-5, check_theta = FALSE, LP_step = NULL, beta = 0.5, gamma = 1, reset_step = FALSE, LP_max_iter = 500, UP_tol = 1e-4, UP_max_iter = 20, startQR = FALSE, verbose = FALSE) {
  if (beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
  if (gamma < 1) stop("Beta must be a nondecreasing factor >= 1")
  if (LP_max_iter < 0 || UP_max_iter < 0) stop("Number of iterations cannot be negative")

  list(method = method, LP_tol_ll = LP_tol_ll, LP_tol_theta = LP_tol_theta, check_theta = check_theta, LP_step = LP_step, beta = beta, gamma = gamma, reset_step = reset_step, LP_max_iter = as.integer(LP_max_iter), UP_tol = UP_tol, UP_max_iter = as.integer(UP_max_iter), startQR = startQR, verbose = verbose)
}
