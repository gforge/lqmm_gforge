#' @title Linear Quantile Models and Linear Quantile Mixed Models
#'
#' @description
#' Fit quantile regression models for independent and hierarchical data
#'
#' \tabular{ll}{ Package: \tab lqmm\cr Type: \tab Package\cr Version: \tab
#' 1.5.5\cr Date: \tab 2019-12-12\cr License: \tab GPL (>=2)\cr LazyLoad: \tab
#' yes\cr }
#'
#' @name lqmm-package
#' @docType package
#' @author Marco Geraci
#'
#' Maintainer: Marco Geraci <geraci@@mailbox.sc.edu>
#' @references Geraci M (2014). Linear quantile mixed models: The lqmm package
#' for Laplace quantile regression. Journal of Statistical Software, 57(13),
#' 1--29. <doi:10.18637/jss.v057.i13>
#'
#' Geraci M and Bottai M (2007). Quantile regression for longitudinal data
#' using the asymmetric Laplace distribution. Biostatistics 8(1), 140--154.
#' <doi:10.1093/biostatistics/kxj039>
#'
#' Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and
#' Computing, 24(3), 461--479. <doi:10.1007/s11222-013-9381-9>.
#' @keywords quantile regression
NULL

#' Labor Pain Data
#'
#' The \code{labor} data frame has 358 rows and 4 columns of the change in pain
#' over time for several 83 women in labor.
#'
#' The labor pain data were reported by Davis (1991) and successively analyzed
#' by Jung (1996) and Geraci and Bottai (2007). The data set consists of
#' repeated measurements of self--reported amount of pain on N = 83 women in
#' labor, of which 43 were randomly assigned to a pain medication group and 40
#' to a placebo group. The response was measured every 30 min on a 100--mm
#' line, where 0 means no pain and 100 means extreme pain. A nearly monotone
#' pattern of missing data was found for the response variable and the maximum
#' number of measurements for each woman was six.
#'
#' @format This data frame contains the following columns: \describe{
#' \item{subject}{ an ordered factor indicating the subject on which the
#' measurement was made.  The levels are labelled \code{1} to \code{83}.  }
#' \item{pain}{ a numeric vector of self--reported pain scores on a 100mm line.
#' } \item{treatment}{ a dummy variable with values \code{1} for subjects who
#' received a pain medication and \code{0} for subjects who received a placebo.
#' } \item{time}{ a numeric vector of times (minutes since randomization) at
#' which \code{pain} was measured.  }
#'
#' }
#' @references Geraci M and Bottai M (2007). Quantile regression for
#' longitudinal data using the asymmetric Laplace distribution. Biostatistics
#' 8(1), 140--154.
#'
#' Jung S (1996). Quasi--likelihood for median regression models. Journal of
#' the American Statistical Association 91, 251--7.
#' @source Davis CS (1991). Semi--parametric and non--parametric methods for
#' the analysis of repeated measurements with applications to clinical trials.
#' Statistics in Medicine 10, 1959--80.
#' @keywords datasets
#' @docType data
#' @name labor
#' @usage data(labor)
NULL

#' Internal lqmm objects
#'
#' Internal lqmm objects.
#'
#' These are not to be called by the user.
#'
#' @aliases .First.lib .Last.lib bandwidth.rq errorHandling extractAll
#' loglik.al loglik.s loglik.t logliki score.al gradientSi switch_check
#' theta.z.dim createLaguerre quad invTfun Tfun F.lqm addnoise permutations
#' asOneFormula allVarsRec nlloglikh rePred C_gradientSh C_gradientSi C_ll_h
#' @keywords internal
#' @name internal
NULL

#' Growth curve data on an orthdontic measurement
#'
#' The \code{Orthodont} data frame has 108 rows and 4 columns of the change in
#' an orthdontic measurement over time for several young subjects.
#'
#' Investigators at the University of North Carolina Dental School followed the
#' growth of 27 children (16 males, 11 females) from age 8 until age 14.  Every
#' two years they measured the distance between the pituitary and the
#' pterygomaxillary fissure, two points that are easily identified on x-ray
#' exposures of the side of the head.
#'
#' @format This data frame contains the following columns: \describe{
#' \item{distance}{ a numeric vector of distances from the pituitary to the
#' pterygomaxillary fissure (mm).  These distances are measured on x-ray images
#' of the skull.  } \item{age}{ a numeric vector of ages of the subject (yr).
#' } \item{Subject}{ an ordered factor indicating the subject on which the
#' measurement was made.  The levels are labelled \code{M01} to \code{M16} for
#' the males and \code{F01} to \code{F13} for the females.  The ordering is by
#' increasing average distance within sex.  } \item{Sex}{ a factor with levels
#' \code{Male} and \code{Female} } }
#' @source Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models
#' in S and S-PLUS}, Springer, New York.  (Appendix A.17)
#'
#' Potthoff, R. F. and Roy, S. N. (1964), ``A generalized multivariate analysis
#' of variance model useful especially for growth curve problems'', Biometrika,
#' 51, 313--326.
#'
#' Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the R
#' Development Core Team (2011). nlme: Linear and Nonlinear Mixed Effects
#' Models. R package version 3.1-100.
#' \url{https://CRAN.R-project.org/package=nlme}
#' @keywords datasets
#' @docType data
#' @name Orthodont
#' @usage data(Orthodont)
NULL
