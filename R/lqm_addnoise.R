##########################################################################################
# QR for counts (Machado, Santos Silva)

addnoise <- function(x, centered = TRUE, B = 0.999) {
  n <- length(x)
  if (centered) {
    z <- x + runif(n, -B / 2, B / 2)
  } else {
    z <- x + runif(n, 0, B)
  }

  return(z)
}
