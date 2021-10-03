lqmmControl <- function(method = "gs", LP_tol_ll = 1e-5, LP_tol_theta = 1e-5, check_theta = FALSE, LP_step = NULL, beta = 0.5, gamma = 1, reset_step = FALSE, LP_max_iter = 500, UP_tol = 1e-4, UP_max_iter = 20, startQR = FALSE, verbose = FALSE) {
  if (beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
  if (gamma < 1) stop("Beta must be a nondecreasing factor >= 1")
  if (LP_max_iter < 0 || UP_max_iter < 0) stop("Number of iterations cannot be negative")

  list(method = method, LP_tol_ll = LP_tol_ll, LP_tol_theta = LP_tol_theta, check_theta = check_theta, LP_step = LP_step, beta = beta, gamma = gamma, reset_step = reset_step, LP_max_iter = as.integer(LP_max_iter), UP_tol = UP_tol, UP_max_iter = as.integer(UP_max_iter), startQR = startQR, verbose = verbose)
}
