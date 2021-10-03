#' Test for Positive Definiteness
#' 
#' This function tests whether all eigenvalues of a symmetric matrix are
#' positive. See \code{help("is.positive.definite")} from package
#' \code{corpcor}.
#' 
#' 
#' @author Original version by Korbinian Strimmer
#' @source Juliane Schaefer, Rainer Opgen-Rhein, Verena Zuber, A. Pedro Duarte
#' Silva and Korbinian Strimmer. (2011). corpcor: Efficient Estimation of
#' Covariance and (Partial) Correlation. R package version 1.6.0.
#' \url{https://CRAN.R-project.org/package=corpcor}
#' @keywords positive definite covariance
is.positive.definite <- function(m, tol, method = c("eigen", "chol")) {
  # source package `corpcor' version 1.6.0 (Juliane Schaefer, Rainer Opgen-Rhein, Verena Zuber, A. Pedro Duarte Silva and Korbinian Strimmer)

  method <- match.arg(method)
  if (!is.matrix(m)) {
    m <- as.matrix(m)
  }
  if (method == "eigen") {
    eval <- eigen(m, only.values = TRUE)$values
    if (is.complex(eval)) {
      warning("Input matrix has complex eigenvalues!")
      return(FALSE)
    }
    if (missing(tol)) {
      tol <- max(dim(m)) * max(abs(eval)) * .Machine$double.eps
    }
    if (sum(eval > tol) == length(eval)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  if (method == "chol") {
    val <- try(chol(m), silent = TRUE)
    if (inherits(val, "try-error")) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}



#' Compute Nearest Positive Definite Matrix
#' 
#' This function computes the nearest positive definite of a real symmetric
#' matrix. See \code{help("make.positive.definite")} from package
#' \code{corpcor}.
#' 
#' 
#' @author Original version by Korbinian Strimmer
#' @source Juliane Schaefer, Rainer Opgen-Rhein, Verena Zuber, A. Pedro Duarte
#' Silva and Korbinian Strimmer. (2011). corpcor: Efficient Estimation of
#' Covariance and (Partial) Correlation. R package version 1.6.0.
#' \url{https://CRAN.R-project.org/package=corpcor}
#' @keywords positive definite covariance
make.positive.definite <- function(m, tol) {
  # source package `corpcor' version 1.6.0 (Juliane Schaefer, Rainer Opgen-Rhein, Verena Zuber, A. Pedro Duarte Silva and Korbinian Strimmer)

  if (!is.matrix(m)) {
    m <- as.matrix(m)
  }
  d <- dim(m)[1]
  if (dim(m)[2] != d) {
    stop("Input matrix is not square!")
  }
  es <- eigen(m)
  esv <- es$values
  if (missing(tol)) {
    tol <- d * max(abs(esv)) * .Machine$double.eps
  }
  delta <- 2 * tol
  tau <- pmax(0, delta - esv)
  dm <- es$vectors %*% diag(tau, d) %*% t(es$vectors)
  return(m + dm)
}
