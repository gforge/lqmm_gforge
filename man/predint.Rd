% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{predint}
\alias{predint}
\title{Predict confidence interval from an `lqmm` object}
\usage{
predint(
  object,
  level = 0,
  alpha = 0.05,
  R = 50,
  seed = round(runif(1, 1, 10000)),
  newdata
)
}
\arguments{
\item{object}{The model fit}

\item{level}{an optional integer vector giving the level of grouping to be
used in obtaining the predictions.}

\item{alpha}{1-\code{alpha} is the confidence level.}

\item{R}{number of bootstrap replications.}

\item{seed}{optional random number generator seed.}

\item{newdata}{Defaults to the original data used for generating the model.
Currently only implemented for `level = 0` predictions.}
}
\value{
A data frame or a list of data frames for \code{predint.lqmm} containing
predictions, lower and upper bounds of prediction intervals, and standard
errors.
}
\description{
A confidence interval predictor for an [lqmm] model fit ([lqm] not implemented)
at the moment.
}
\examples{

## Orthodont data
data(Orthodont)

# Random intercept model
fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	tau = c(0.1,0.5,0.9), data = Orthodont)

# 95\% confidence intervals
predint(fitOi.lqmm, level = 0, alpha = 0.05)
}
\references{
Geraci M and Bottai M (2014). Linear quantile mixed models.
Statistics and Computing, 24(3), 461--479.
}
\seealso{
\code{\link{lqmm}}, \code{\link{ranef.lqmm}},
\code{\link{coef.lqmm}}
}
\author{
Marco Geraci
}
\keyword{prediction}
