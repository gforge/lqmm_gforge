residuals.lqmm <- function(object, level = 0, ...) {
  object$y[object$revOrder] - predict(object, level = level)
}
