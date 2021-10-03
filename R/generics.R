##########################################################################################
# New generics



#' Bootstrap functions for LQM and LQMM
#' 
#' This function is used to obtain a bootstrap sample of a fitted LQM or LQMM.
#' It is a generic function.
#' 
#' 
#' @aliases boot.lqm boot.lqmm boot
#' @param object an object of \code{\link{class}} "lqm" or "lqmm".
#' @param R number of bootstrap replications.
#' @param seed optional random number generator seed.
#' @param startQR logical flag. If \code{TRUE} the estimated parameters in
#' \code{object} are used as starting values in the algorithm applied to each
#' bootstrap sample. This may cause the algorithm to converge too often to a
#' similar optimum, which would ultimately result in underestimated standard
#' errors. If \code{FALSE} (recommended), starting values are based on
#' \code{\link{lm}}.
#' @return An object of class \code{boot.lqm} is a data frame with \code{R}
#' rows and \code{npars} columns containing the bootstrap estimates of
#' \code{theta}. If \code{object} contains results for multiple quantiles,
#' \code{boot.lqm} returns an array of dimension \code{c(R,npars,nt)}, where
#' \code{nt} is the length of \code{tau}.
#' 
#' An object of class \code{boot.lqmm} is a data frame with \code{R} rows and
#' \code{npars} columns containing the bootstrap estimates of \code{theta_x},
#' \code{theta_z}, and \code{scale}. If \code{object} contains results for
#' multiple quantiles, \code{boot.lqmm} returns an array of dimension
#' \code{c(R,npars,nt)}, where \code{nt} is the length of \code{tau}. The
#' elements of \code{theta_z} are labelled with \code{reStruct}. See function
#' \code{\link{covHandling}} and the example below on how to derive the
#' variance-covariance matrix of the random effects starting from
#' \code{theta_z}.
#' 
#' The following attributes are available:
#' 
#' \item{tau}{index of the quantile(s).} \item{estimated}{the estimated
#' parameter as given by \code{object}.} \item{R}{number of bootstrap
#' replications.} \item{seed}{the random number generator seed used to produce
#' the bootstrap sample.} \item{npars}{total numer of parameters.}
#' \item{rdf}{the number of residual degrees of freedom.} \item{indices}{the
#' bootstrap sample of independent data units.}
#' @author Marco Geraci
#' @keywords bootstrap standard errors
#' @examples
#' 
#' 
#' # boot.lqm
#' set.seed(123)
#' n <- 500
#' test <- data.frame(x = runif(n,0,1))
#' test$y <- 30 + test$x + rnorm(n)
#' fit.lqm <- lqm(y ~ x, data = test, tau = 0.5)
#' fit.boot <- boot(fit.lqm)
#' str(fit.boot)
#' 
#' # boot.lqmm
#' data(Orthodont)
#' fit <- lqmm(distance ~ age, random = ~ 1, group = Subject,
#' 	tau = 0.5, data = Orthodont)
#' fit.boot <- boot(fit)
#' str(fit.boot)
#' 
boot <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE) UseMethod("boot")


#' Extract Fixed and Random Bootstrapped Parameters
#' 
#' This generic function extracts the fixed and random components of
#' bootstrapped estimates of an \code{lqmm} object.
#' 
#' The \code{"random"} parameters refer to the "raw" parameters of the
#' variance-covariance matrix of the random effects as returned by
#' \code{\link{lqmm.fit.gs}} and \code{\link{lqmm.fit.df}}.
#' 
#' @aliases extractBoot.boot.lqmm extractBoot
#' @param object an object of \code{\link{class}} \code{boot.lqmm}.
#' @param which character indicating whether \code{"fixed"} or \code{"random"}
#' parameters.
#' @return a matrix of bootstrapped estimates.
#' @author Marco Geraci
#' @seealso \code{\link{boot.lqmm}}, \code{\link{lqmm.fit.gs}},
#' \code{\link{lqmm.fit.df}}
#' @keywords bootstrap
#' @examples
#' 
#' 
#' ## Orthodont data
#' data(Orthodont)
#' 
#' # Random intercept model
#' fit <- lqmm(distance ~ age, random = ~ 1, group = Subject,
#' 	tau = 0.5, data = Orthodont)
#' fit.boot <- boot(fit)
#' 
#' # extract fixed effects
#' B <- extractBoot(fit.boot, which = "fixed")
#' 
#' # covariance matrix estimated fixed parameters
#' cov(B)
#' 
#' 
extractBoot <- function(object, which = "fixed") UseMethod("extractBoot")
predint <- function(object, level = 0, alpha = 0.05, R = 50, seed = round(runif(1, 1, 10000))) UseMethod("predint")
