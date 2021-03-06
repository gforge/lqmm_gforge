% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lqm_S3_coef.R
\name{coef.lqm}
\alias{coef.lqm}
\alias{coef.lqm.counts}
\title{Extract LQM Coefficients}
\usage{
\method{coef}{lqm}(object, ...)
}
\arguments{
\item{object}{an \code{lqm} or \code{lqm.counts} object.}

\item{\dots}{not used.}
}
\value{
a vector for single quantiles or a matrix for multiple quantiles.
}
\description{
\code{coef} extracts model coefficients from \code{lqm}, \code{lqm.counts}
objects.
}
\seealso{
\code{\link{lqm}} \code{\link{summary.lqm}}
\code{\link{lqm.counts}}
}
\author{
Marco Geraci
}
\keyword{coefficients}
