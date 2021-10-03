residuals.lqm <- function(object, ...) {
  ans <- as.numeric(object$y) - predict(object)
  return(ans)
}
