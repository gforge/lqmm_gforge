#' Extract Variance-Covariance Matrix
#' 
#' This function extracts the variance-covariance matrix of the random effects
#' from a fitted \code{lqmm} object.
#' 
#' This function returns the variance or the variance-covariance matrix of the
#' random effects. It calls \code{\link{covHandling}} to manage the output of
#' \code{\link{lqmm.fit.gs}} or \code{\link{lqmm.fit.df}}. A post-fitting
#' approximation to the nearest positive (semi)definite matrix (Higham, 2002)
#' is applied if necessary. The generic function \code{VarCorr} is imported
#' from the \code{nlme} package (Pinheiro et al, 2014).
#' 
#' @aliases VarCorr VarCorr.lqmm
#' @param x an object of \code{\link{class}} "lqmm".
#' @param sigma not used.
#' @param ...  not used.
#' @author Marco Geraci
#' @seealso \code{\link{lqmm}} \code{\link{coef.lqmm}}
#' @references Higham N (2002). Computing the Nearest Correlation Matrix - A
#' Problem from Finance. IMA Journal of Numerical Analysis, 22, 329-343.
#' 
#' Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear
#' and Nonlinear Mixed Effects Models. R package version 3.1-117,
#' \url{https://CRAN.R-project.org/package=nlme}.
#' @keywords covariance coefficients
VarCorr.lqmm <- function(x, sigma = NULL, ...) {
  tau <- x$tau
  nq <- length(tau)

  theta_z <- x$theta_z
  dim_theta <- x$dim_theta
  q <- dim_theta[2]
  cov_name <- x$cov_name
  type <- x$type
  mm <- x$mm

  if (nq == 1) {
    sigma <- covHandling(theta = theta_z, n = q, cov_name = cov_name, quad_type = type)
    if (cov_name == "pdIdent") {
      sigma <- rep(sigma, q)
      names(sigma) <- mm
    }
    if (cov_name == "pdDiag") {
      names(sigma) <- mm
    }
    if (cov_name == "pdCompSymm") {
      rownames(sigma) <- colnames(sigma) <- mm
    }
    if (cov_name == "pdSymm") {
      rownames(sigma) <- colnames(sigma) <- mm
    }
  } else {
    sigma <- vector("list", nq)
    names(sigma) <- format(tau, digits = 4)
    for (i in 1:nq) {
      sigma[[i]] <- covHandling(theta = theta_z[, i], n = q, cov_name = cov_name, quad_type = type)
      if (cov_name == "pdIdent") {
        sigma[[i]] <- rep(sigma[[i]], q)
        names(sigma[[i]]) <- mm
      }
      if (cov_name == "pdDiag") {
        names(sigma[[i]]) <- mm
      }
      if (cov_name == "pdCompSymm") {
        rownames(sigma[[i]]) <- colnames(sigma[[i]]) <- mm
      }
      if (cov_name == "pdSymm") {
        rownames(sigma[[i]]) <- colnames(sigma[[i]]) <- mm
      }
    }
  }

  return(sigma)
}
