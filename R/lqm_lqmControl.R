#' Control parameters for lqm estimation
#' 
#' A list of parameters for controlling the fitting process.
#' 
#' The methods "gs1" and "gs2" implement the same algorithm (Bottai et al,
#' 2015). The former is based on C code, the latter on R code. While the C code
#' is faster, the R code seems to be more efficient in handling large datasets.
#' For method "gs2", it is possible to replace the classical non-differentiable
#' loss function with a smooth version (Chen, 2007).
#' 
#' @param method character vector that specifies which code to use for carrying
#' out the gradient search algorithm: "gs1" (default) based on C code and "gs2"
#' based on R code. Method "gs3" uses a smoothed loss function. See details.
#' @param loop_tol_ll tolerance expressed as relative change of the
#' log-likelihood.
#' @param loop_tol_theta tolerance expressed as relative change of the
#' estimates.
#' @param check_theta logical flag. If \code{TRUE} the algorithm performs a
#' check on the change in the estimates in addition to the likelihood.
#' @param loop_step step size (default standard deviation of response).
#' @param beta decreasing step factor for line search (0,1).
#' @param gamma nondecreasing step factor for line search (>= 1).
#' @param reset_step logical flag. If \code{TRUE} the step size is re-setted to
#' the initial value at each iteration.
#' @param loop_max_iter maximum number of iterations.
#' @param smooth logical flag. If \code{TRUE} the standard loss function is
#' replaced with a smooth approximation.
#' @param omicron small constant for smoothing the loss function when using
#' \code{smooth = TRUE}. See details.
#' @param verbose logical flag.
#' @return a list of control parameters.
#' @author Marco Geraci
#' @seealso \code{\link{lqm}}
#' @references Bottai M, Orsini N, Geraci M (2015). A Gradient Search
#' Maximization Algorithm for the Asymmetric Laplace Likelihood, Journal of
#' Statistical Computation and Simulation, 85(10), 1919-1925.
#' 
#' Chen C (2007). A finite smoothing algorithm for quantile regression. Journal
#' of Computational and Graphical Statistics, 16(1), 136-164.
#' @keywords fitting control
lqmControl <- function(method = "gs1", loop_tol_ll = 1e-5, loop_tol_theta = 1e-3, check_theta = FALSE, loop_step = NULL, beta = 0.5, gamma = 1.25, reset_step = FALSE, loop_max_iter = 1000, smooth = FALSE, omicron = 0.001, verbose = FALSE) {
  if (beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
  if (gamma < 1) stop("Beta must be a nondecreasing factor >= 1")
  if (loop_max_iter < 0) stop("Number of iterations cannot be negative")

  list(method = method, loop_tol_ll = loop_tol_ll, loop_tol_theta = loop_tol_theta, check_theta = check_theta, loop_step = loop_step, beta = beta, gamma = gamma, reset_step = reset_step, loop_max_iter = as.integer(loop_max_iter), smooth = smooth, omicron = omicron, verbose = verbose)
}
