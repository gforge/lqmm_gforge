#' Fitting Linear Quantile Mixed Models
#'
#' \code{lqmm} is used to fit linear quantile mixed models based on the
#' asymmetric Laplace distribution.
#'
#' The function computes an estimate on the tau-th quantile function of the
#' response, conditional on the covariates, as specified by the \code{formula}
#' argument, and on random effects, as specified by the \code{random} argument.
#' The quantile predictor is assumed to be linear. The function maximizes the
#' (log)likelihood of the Laplace regression proposed by Geraci and Bottai
#' (2014). The likelihood is numerically integrated via Gaussian quadrature
#' techniques. The optimization algorithm is based on the gradient of the
#' Laplace log--likelihood (\code{control = list(method = "gs")}). An
#' alternative optimization algorithm is based on a Nelder-Mead algorithm
#' (\code{control = list(method = "df")}) via \code{\link{optim}}. The scale
#' parameter is optimized in a refinement step via \code{\link{optimize}}.
#'
#' Quadrature approaches include Gauss-Hermite (\code{type = "normal"}) and
#' Gauss-Laguerre (\code{type = "robust"}) quadrature. The argument \code{rule}
#' takes one of the following: 1 (product rule quadrature), 2 (sparse grid
#' quadrature), 3 (nested quadrature rule - only for \code{type = "normal"}), 4
#' (quadrature rule with the smallest number of nodes between rules 1 or 2).
#' Rules 2 and 3 have not yet been tested extensively.
#'
#' Different standard types of positive--definite matrices for the random
#' effects can be specified: \code{pdIdent} multiple of an identity;
#' \code{pdCompSymm} compound symmetry structure (constant diagonal and
#' constant off--diagonal elements); \code{pdDiag} diagonal; \code{pdSymm}
#' general positive--definite matrix, with no additional structure.
#'
#' Weights are given to clusters, therefore it is expected that these are
#' constant within cluster. When the weights are specified in the main call,
#' then the first value by \code{group} in the vector \code{weights} will be
#' replicated for the same length of each group. Alternatively, different
#' weights within the same cluster can be introduced with a direct call to
#' \code{\link{lqmm.fit.gs} or \link{lqmm.fit.df}}.
#'
#' The \code{lqmm} vignette can be accessed by typing \code{help(package =
#' "lqmm")} and then following the link 'User guides, package vignettes and
#' other documentation'.
#'
#' @param fixed an object of class \code{\link{formula}} for fixed effects: a
#' symbolic description of the model to be fitted.
#' @param random a one-sided formula of the form \code{~x1 + x2 + ... + xn} for
#' random effects: a symbolic description of the model to be fitted.
#' @param group grouping factor.
#' @param covariance variance--covariance matrix of the random effects. Default
#' is \code{pdDiag} (see details).
#' @param tau the quantile(s) to be estimated.
#' @param nK number of quadrature knots.
#' @param type type of quadrature "c("normal","robust")" (see details).
#' @param rule quadrature rule (see details).
#' @param data an optional data frame containing the variables named in
#' \code{fixed}, \code{random} and \code{group}. By default the variables are
#' taken from the environment from which \code{lqmm} is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#' process of the same length as the number of rows of \code{data}. Weights are
#' given to clusters, therefore units within the same cluster receive the same
#' weight (see details).
#' @param na.action a function that indicates what should happen when the data
#' contain \code{NA}s.  The default action (\code{na.fail}) causes \code{lqmm}
#' to print an error message and terminate if there are any incomplete
#' observations.
#' @param control list of control parameters of the fitting process. See
#' \code{\link{lqmmControl}}.
#' @param contrasts not yet implemented.
#' @param fit logical flag. If FALSE the function returns a list of arguments
#' to be passed to \code{lqmm.fit}.
#' @return \code{lqmm} returns an object of \code{\link{class}} \code{lqmm}.
#'
#' The function \code{summary} is used to obtain and print a summary of the
#' results.
#'
#' An object of class \code{lqmm} is a list containing the following
#' components:
#'
#' \item{theta}{a vector containing fixed regression coefficients and
#' parameters of the variance--covariance matrix of the random effects. See
#' \code{\link{VarCorr.lqmm}} to extract the variance--covariance of the random
#' effects from an "lqmm" object.} \item{theta_x,theta_z}{partition of
#' \code{theta}: fixed regression coefficients (\code{theta_x}) and unique
#' variance--covariance parameters (\code{theta_z}).} \item{scale}{the scale
#' parameter.} \item{gradient}{the gradient (\code{control = list(method =
#' "gs")}).} \item{logLik}{the log--likelihood.} \item{opt}{details on
#' optimization (see \code{\link{lqmm.fit.gs}} and \code{\link{lqmm.fit.df}}).}
#' \item{call}{the matched call.} \item{nn}{column names of \code{mmf}.}
#' \item{mm}{column names of \code{mmr}.} \item{nobs}{the number of
#' observations.} \item{dim_theta}{the number of columns in \code{mmf} and
#' \code{mmr}.} \item{dim_theta_z}{the length of \code{theta_z}.}
#' \item{edf}{length of \code{theta}.} \item{rdf}{the number of residual
#' degrees of freedom.} \item{df}{edf + 1 (scale parameter).} \item{tau}{the
#' estimated quantile(s).} \item{mmf}{the model matrix -- fixed effects.}
#' \item{mmr}{the model matrix -- random effects.} \item{y}{the model
#' response.} \item{revOrder}{original order of observations (now ordered
#' according to \code{group}).} \item{weights}{the likelihood weights used in
#' the fitting process (a vector of 1's if \code{weights} is missing or
#' \code{NULL}).} \item{group}{the grouping factor.} \item{ngroups}{the number
#' of groups.} \item{QUAD}{quadrature nodes and weights.} \item{type}{the type
#' of quadrature.} \item{rule}{quadrature rule.} \item{InitialPar}{starting
#' values for theta.} \item{control}{list of control parameters used for
#' optimization (see \code{\link{lqmmControl}}).} \item{cov_name}{class of
#' variance-covariance matrix for the random effects.} \item{mfArgs}{arguments
#' for \code{\link{model.frame}} to return the full data frame.}
#' @note Updates/FAQ/news are published here
#' \url{http://marcogeraci.wordpress.com/}. New versions are usually published
#' here \url{https://r-forge.r-project.org/R/?group_id=1396} before going on
#' CRAN.
#' @author Marco Geraci
#' @seealso \code{\link{lqm}, \link{summary.lqmm}, \link{coef.lqmm},
#' \link{VarCorr.lqmm}, \link{predict.lqmm}, \link{residuals.lqmm}}
#' @references Genz A, and Keister BD (1996). Fully symmetric interpolatory
#' rules for multiple integrals over infinite regions with Gaussian weight.
#' Journal of Computational and Applied Mathematics, 71(2), 299--309.
#' <doi:10.1016/0377-0427(95)00232-4>
#'
#' Geraci M (2014). Linear quantile mixed models: The lqmm package for Laplace
#' quantile regression. Journal of Statistical Software, 57(13), 1--29.
#' <doi:10.18637/jss.v057.i13>
#'
#' Geraci M and Bottai M (2007). Quantile regression for longitudinal data
#' using the asymmetric Laplace distribution. Biostatistics 8(1), 140--154.
#' <doi:10.1093/biostatistics/kxj039>
#'
#' Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and
#' Computing, 24(3), 461--479. <doi:10.1007/s11222-013-9381-9>.
#'
#' Heiss F, and Winschel V (2008). Likelihood approximation by numerical
#' integration on sparse grids. Journal of Econometrics, 144(1), 62--80.
#' <doi:10.1016/j.jeconom.2007.12.004>
#' @keywords quantile regression
#' @examples
#'
#'
#' # Test example
#' set.seed(123)
#'
#' M <- 50
#' n <- 10
#' test <- data.frame(x = runif(n*M,0,1), group = rep(1:M,each=n))
#' test$y <- 10*test$x + rep(rnorm(M, 0, 2), each = n) + rchisq(n*M, 3)
#' fit.lqmm <- lqmm(fixed = y ~ x, random = ~ 1, group = group,
#' 	data = test, tau = 0.5, nK = 11, type = "normal")
#' fit.lqmm
#'
#' #Call: lqmm(fixed = y ~ x, random = ~1, group = group, tau = 0.5, nK = 11,
#' #    type = "normal", data = test)
#' #Quantile 0.5
#'
#' #Fixed effects:
#' #(Intercept)            x
#' #      3.443        9.258
#'
#' #Covariance matrix of the random effects:
#' #(Intercept)
#' #      3.426
#'
#' #Residual scale parameter: 0.8697 (standard deviation 2.46)
#' #Log-likelihood: -1178
#'
#' #Number of observations: 500
#' #Number of groups: 50
#'
#'
#' ## Orthodont data
#' data(Orthodont)
#'
#' # Random intercept model
#' fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
#' 	tau = c(0.1,0.5,0.9), data = Orthodont)
#' coef(fitOi.lqmm)
#'
#' # Random slope model
#' fitOs.lqmm <- lqmm(distance ~ age, random = ~ age, group = Subject,
#' 	tau = c(0.1,0.5,0.9), cov = "pdDiag", data = Orthodont)
#'
#' # Extract estimates
#' VarCorr(fitOs.lqmm)
#' coef(fitOs.lqmm)
#' ranef(fitOs.lqmm)
#'
#' # AIC
#' AIC(fitOi.lqmm)
#' AIC(fitOs.lqmm)
lqmm <- function(fixed, random, group, covariance = "pdDiag", tau = 0.5, nK = 7, type = "normal", rule = 1, data = sys.frame(sys.parent()), subset, weights, na.action = na.fail, control = list(), contrasts = NULL, fit = TRUE) {
  Call <- match.call()

  if (any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
  nq <- length(tau)

  if (!is.data.frame(data)) stop("`data' must be a data frame")
  if (!(type %in% c("normal", "robust"))) stop("type must be either `normal' or `robust'")


  # check arguments
  if (!inherits(fixed, "formula") || length(fixed) != 3) {
    stop("\nFixed-effects model must be a formula of the form \"y ~ x\"")
  }
  if (!inherits(random, "formula") || length(random) != 2) {
    stop("\nRandom-effects model must be a formula of the form \" ~ x\"")
  }

  groupFormula <- asOneSidedFormula(Call[["group"]])
  group <- groupFormula[[2]]

  ## extract data frame with all necessary information

  mfArgs <- list(formula = asOneFormula(random, fixed, group), data = data, na.action = na.action)
  if (!missing(subset)) {
    mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
  }
  if (!missing(weights)) {
    mfArgs[["weights"]] <- weights
  }
  mfArgs$drop.unused.levels <- TRUE
  dataMix <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMix) # preserve the original order
  for (i in names(contrasts)) contrasts(dataMix[[i]]) <- contrasts[[i]]

  ## sort the model.frame by groups
  grp <- model.frame(groupFormula, dataMix)

  ## ordering data by groups
  ord <- order(unlist(grp, use.names = FALSE))
  grp <- grp[ord, , drop = TRUE]
  dataMix <- dataMix[ord, , drop = FALSE]
  revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
  ngroups <- length(unique(grp))

  ## obtaining model matrices and vectors
  y <- eval(fixed[[2]], dataMix)
  mmr <- model.frame(random, dataMix)
  mmr <- model.matrix(random, data = mmr)
  # Likelihood weights
  if (!missing(weights)) weights <- model.weights(dataMix)[!duplicated(grp)]
  if (!missing(weights) && is.null(weights)) weights <- rep(1, ngroups)
  if (missing(weights)) weights <- rep(1, ngroups)
  # keeping the contrasts for use in predict
  contr <- attr(mmr, "contr")
  mmf <- model.frame(fixed, dataMix)
  Terms <- attr(mmf, "terms")
  auxContr <- lapply(mmf, function(el) {
    if (inherits(el, "factor") && length(levels(el)) > 1) contrasts(el)
  })

  contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
  contr <- contr[!unlist(lapply(contr, is.null))]
  mmf <- model.matrix(fixed, data = mmf)

  ## define dimensions
  cov_name <- covariance
  if (type == "robust" & !(cov_name %in% c("pdIdent", "pdDiag"))) stop("Gauss-Laguerre quadrature available only for uncorrelated random effects.")

  dim_theta <- integer(2)
  dim_theta[1] <- ncol(mmf)
  dim_theta[2] <- ncol(mmr)
  dim_theta_z <- theta.z.dim(type = cov_name, n = dim_theta[2])

  ## Check if product rule quadrature is computationally heavy

  if (rule == 1) {
    if (dim_theta[2] > 4 && nK > 11) {
      warning(paste("For current value of \"nK\" the total number of quadrature knots is ", nK^dim_theta[2], sep = ""))
    }
  }

  ## Quandrature nodes and weights
  QUAD <- quad(q = dim_theta[2], k = nK, type = type, rule = rule)

  ## Control
  if (is.null(names(control))) {
    control <- lqmmControl()
  } else {
    control_default <- lqmmControl()
    control_names <- intersect(names(control), names(control_default))
    control_default[control_names] <- control[control_names]
    control <- control_default
  }
  if (is.null(control$LP_step)) control$LP_step <- sd(as.numeric(y))
  method <- control$method
  if (method == "gs" & cov_name == "pdCompSymm") {
    method <- "df"
    cat("Switching to Nelder-Mead optimization \n")
  }

  # Starting values
  theta_z <- if (type == "normal") rep(1, dim_theta_z) else rep(invvarAL(1, 0.5), dim_theta_z)
  lmfit <- lm.wfit(x = mmf, y = y, w = rep(weights, table(grp)))
  theta_x <- lmfit$coefficients

  if (control$startQR) {
    q_close <- if (nq == 1) tau else 0.5
    fit_rq <- lqm.fit.gs(
      theta = theta_x, x = as.matrix(mmf), y = y, weights = rep(weights, table(grp)), tau = q_close,
      control = lqmControl(loop_step = sd(as.numeric(y)))
    )
    theta_x <- fit_rq$theta
    sigma_0 <- fit_rq$scale
  } else {
    sigma_0 <- invvarAL(mean(lmfit$residuals^2), 0.5)
  }
  theta_0 <- c(theta_x, theta_z)

  ## Create list with all necessary arguments
  FIT_ARGS <- list(theta_0 = theta_0, x = as.matrix(mmf), y = y, z = as.matrix(mmr), weights = weights, V = QUAD$nodes, W = QUAD$weights, sigma_0 = sigma_0, tau = tau, group = grp, cov_name = cov_name, control = control)

  if (!fit) {
    return(FIT_ARGS)
  }

  ## Estimation
  if (method == "gs") {
    if (nq == 1) {
      fit <- do.call(lqmm.fit.gs, FIT_ARGS)
    } else {
      fit <- vector("list", nq)
      names(fit) <- format(tau, digits = 4)
      for (i in 1:nq) {
        FIT_ARGS$tau <- tau[i]
        fit[[i]] <- do.call(lqmm.fit.gs, FIT_ARGS)
      }
    }
  }
  if (method == "df") {
    if (nq == 1) {
      fit <- do.call(lqmm.fit.df, FIT_ARGS)
    } else {
      fit <- vector("list", nq)
      names(fit) <- format(tau, digits = 4)
      for (i in 1:nq) {
        FIT_ARGS$tau <- tau[i]
        fit[[i]] <- do.call(lqmm.fit.df, FIT_ARGS)
      }
    }
  }

  nn <- colnames(mmf)
  mm <- colnames(mmr)

  if (nq > 1) {
    fit$theta_x <- matrix(NA, dim_theta[1], nq)
    fit$theta_z <- matrix(NA, dim_theta_z, nq)
    for (i in 1:nq) {
      fit$theta_x[, i] <- fit[[i]]$theta_x <- fit[[i]]$theta[1:dim_theta[1]]
      fit$theta_z[, i] <- fit[[i]]$theta_z <- fit[[i]]$theta[-(1:dim_theta[1])]
    }
    rownames(fit$theta_x) <- nn
    colnames(fit$theta_x) <- colnames(fit$theta_z) <- format(tau, digits = 4)
  } else {
    fit$theta_x <- fit$theta[1:dim_theta[1]]
    fit$theta_z <- fit$theta[-(1:dim_theta[1])]
  }


  fit$call <- Call
  fit$nn <- nn
  fit$mm <- mm
  fit$nobs <- length(y)
  fit$dim_theta <- dim_theta
  fit$dim_theta_z <- dim_theta_z
  fit$edf <- fit$dim_theta[1] + fit$dim_theta_z
  fit$rdf <- fit$nobs - fit$edf
  fit$df <- dim_theta[1] + dim_theta_z + 1
  fit$tau <- tau
  fit$mmf <- as.matrix(mmf)
  fit$mmr <- as.matrix(mmr)
  fit$y <- y
  fit$revOrder <- revOrder
  fit$weights <- weights
  fit$contrasts <- contr
  fit$group <- grp
  fit$ngroups <- ngroups
  fit$QUAD <- QUAD
  fit$type <- type
  fit$rule <- rule
  fit$InitialPar <- list(theta = theta_0, sigma = sigma_0)
  fit$control <- control
  fit$cov_name <- cov_name
  # attr(fit$cov_name, "cov_type") <- cov.sel(cov_name)
  fit$mfArgs <- mfArgs

  class(fit) <- "lqmm"
  fit
}
