##########################################################################################
# QR for counts (Machado, Santos Silva)

lqm.counts <- function(formula, data, weights = NULL, offset = NULL, contrasts = NULL, tau = 0.5, M = 50, zeta = 1e-5, B = 0.999, cn = NULL, alpha = 0.05, control = list()) {
  nq <- length(tau)
  if (nq > 1) {
    stop("One quantile at a time")
  }

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) {
    stop("'weights' must be a numeric vector")
  }
  if (is.null(w)) {
    w <- rep(1, length(y))
  }
  x <- model.matrix(mt, mf, contrasts)
  p <- ncol(x)
  n <- nrow(x)
  term.labels <- colnames(x)

  if (is.null(offset)) {
    offset <- rep(0, n)
  }
  if (is.null(names(control))) {
    control <- lqmControl()
  } else {
    control_default <- lqmControl()
    control_names <- intersect(names(control), names(control_default))
    control_default[control_names] <- control[control_names]
    control <- control_default
  }
  if (is.null(control$loop_step)) {
    control$loop_step <- sd(as.numeric(y))
  }
  if (control$beta > 1 || control$beta < 0) {
    stop("Beta must be a decreasing factor in (0,1)")
  }
  if (control$gamma < 1) {
    stop("Beta must be a nondecreasing factor >= 1")
  }
  if (p == 1) {
    control$loop_tol_ll <- 0.005
  }
  theta_0 <- glm.fit(x = as.matrix(x), y = y, weights = w, offset = offset, family = poisson())$coefficients

  # Add noise
  Z <- replicate(M, addnoise(y, centered = FALSE, B = B))
  # Transform Z
  TZ <- apply(Z, 2, function(x, off, tau, zeta) {
    log(ifelse((x -
      tau) > zeta, x - tau, zeta)) - off
  }, off = offset, tau = tau, zeta = zeta)
  # Fit linear QR on TZ
  fit <- apply(TZ, 2, function(y, x, weights, tau, control,
                               theta) {
    lqm.fit.gs(theta = theta, x = x, y = y, weights = weights, tau = tau, control = control)
  }, x = x, weights = w, tau = tau, control = control, theta = theta_0)
  # Trasform back
  yhat <- sapply(fit, function(obj, x) x %*% obj$theta, x = x)
  yhat <- as.matrix(yhat)
  eta <- sweep(yhat, 1, offset, "+")
  zhat <- tau + exp(eta)
  #
  Fvec <- Vectorize(F.lqm)
  if (is.null(cn)) cn <- 0.5 * log(log(n)) / sqrt(n)
  F <- apply(zhat, 2, Fvec, cn = cn)
  Fp <- apply(zhat + 1, 2, Fvec, cn = cn)

  multiplier <- (tau - (TZ <= yhat))^2
  a <- array(NA, dim = c(p, p, M))
  for (i in 1:M) a[, , i] <- t(x * multiplier[, i]) %*% x / n

  multiplier <- tau^2 + (1 - 2 * tau) * (y <= (zhat - 1)) +
    ((zhat - y) * (zhat - 1 < y & y <= zhat)) * (zhat - y -
      2 * tau)
  b <- array(NA, dim = c(p, p, M))
  for (i in 1:M) b[, , i] <- t(x * multiplier[, i]) %*% x / n

  multiplier <- exp(eta) * (F <= Z & Z < Fp)
  d <- array(NA, dim = c(p, p, M))
  sel <- rep(TRUE, M)
  for (i in 1:M) {
    tmpInv <- try(solve(t(x * multiplier[, i]) %*% x / n), silent = TRUE)
    if (!inherits(tmpInv, "try-error")) {
      d[, , i] <- tmpInv
    } else {
      sel[i] <- FALSE
    }
  }

  dad <- 0
  dbd <- 0
  for (i in (1:M)[sel]) {
    dad <- dad + d[, , i] %*% a[, , i] %*% d[, , i]
    dbd <- dbd + d[, , i] %*% b[, , i] %*% d[, , i]
  }

  m.n <- sum(sel)
  if (m.n != 0) {
    V <- dad / (m.n^2) + (1 - 1 / m.n) * dbd * 1 / m.n
    V <- V / n
    stds <- sqrt(diag(V))
  } else {
    stds <- NA
    warning("Standard error not available")
  }

  est <- sapply(fit, function(x) x$theta)
  est <- if (p == 1) mean(est) else rowMeans(est)

  qfit <- if (p == 1) {
    tau + exp(mean(eta[1, ]))
  } else {
    tau + exp(rowMeans(eta))
  }

  lower <- est + qt(alpha / 2, n - p) * stds
  upper <- est + qt(1 - alpha / 2, n - p) * stds
  tP <- 2 * pt(-abs(est / stds), n - p)

  ans <- cbind(est, stds, lower, upper, tP)
  colnames(ans) <- c(
    "Value", "Std. Error", "lower bound", "upper bound",
    "Pr(>|t|)"
  )
  rownames(ans) <- names(est) <- term.labels

  fit <- list()
  fit$call <- call
  fit$na.action <- attr(mf, "na.action")
  fit$contrasts <- attr(x, "contrasts")
  fit$term.labels <- term.labels
  fit$terms <- mt

  fit$theta <- est
  fit$tau <- tau
  fit$nobs <- n
  fit$M <- M
  fit$Mn <- m.n
  fit$rdf <- n - p
  fit$x <- x
  fit$y <- y
  fit$fitted <- qfit
  fit$offset <- offset
  fit$Cov <- V
  fit$tTable <- ans
  fit$levels <- .getXlevels(mt, mf)
  fit$InitialPar <- list(theta = theta_0)
  fit$control <- control

  class(fit) <- "lqm.counts"

  return(fit)
}
