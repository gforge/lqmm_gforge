% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lqmm_S3_summary.R
\name{summary.lqmm}
\alias{summary.lqmm}
\title{Summary for an \code{lqmm} Object}
\usage{
\method{summary}{lqmm}(object, method = "boot", alpha = 0.05, covariance = FALSE, ...)
}
\arguments{
\item{object}{an object of \code{\link{class}} \code{lqmm}.}

\item{method}{specifies the method used to compute standard errors.
Currently, only the bootstrap method ("boot") is available.}

\item{alpha}{significance level.}

\item{covariance}{logical flag. If \code{TRUE} the covariance matrix of the
bootstrap estimates is provided.}

\item{\dots}{see \code{\link{boot.lqmm}} for additional arguments.}
}
\value{
an object of class \code{summary.lqmm}. The function
\code{summary.lqmm} computes and returns a list of summary statistics of the
fitted linear quantile mixed model given in \code{object}, using the
components (list elements) from its argument, plus

\item{Cov}{the covariance matrix obtained from the bootstrapped estimates
(if \code{covariance = TRUE}).} \item{tTable}{a matrix with estimates,
standard errors, etc.} \item{B}{the matrix of all bootstrapped parameters.}
}
\description{
Summary method for class \code{lqmm}.
}
\details{
\code{print.summary.lqmm} formats the coefficients, standard errors, etc.
and additionally gives `significance stars'.
}
\examples{

data(Orthodont)
fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	tau = c(0.1,0.5,0.9), data = Orthodont)
summary(fitOi.lqmm)
}
\seealso{
\code{\link{print.summary.lqmm}} \code{\link{lqmm}}
}
\author{
Marco Geraci
}
\keyword{bootstrap}
\keyword{errors}
\keyword{standard}
