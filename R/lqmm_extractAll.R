extractAll <- function(object) {
  nq <- length(object$tau)

  dim_theta_z <- object$dim_theta_z
  dimn <- c(object$nn, paste("reStruct", 1:dim_theta_z, sep = ""), "scale")

  if (nq == 1) {
    ans <- c(object$theta, object$scale)
    names(ans) <- dimn
  } else {
    val <- NULL
    for (i in 1:nq) val <- c(val, object[[i]]$scale)
    ans <- rbind(object$theta_x, object$theta_z, val)
    rownames(ans) <- dimn
  }
  return(ans)
}
