lqmControl <- function(method = "gs1", loop_tol_ll = 1e-5, loop_tol_theta = 1e-3, check_theta = FALSE, loop_step = NULL, beta = 0.5, gamma = 1.25, reset_step = FALSE, loop_max_iter = 1000, smooth = FALSE, omicron = 0.001, verbose = FALSE) {
  if (beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
  if (gamma < 1) stop("Beta must be a nondecreasing factor >= 1")
  if (loop_max_iter < 0) stop("Number of iterations cannot be negative")

  list(method = method, loop_tol_ll = loop_tol_ll, loop_tol_theta = loop_tol_theta, check_theta = check_theta, loop_step = loop_step, beta = beta, gamma = gamma, reset_step = reset_step, loop_max_iter = as.integer(loop_max_iter), smooth = smooth, omicron = omicron, verbose = verbose)
}
