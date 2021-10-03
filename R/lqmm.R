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
