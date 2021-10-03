F.lqm <- function(x, cn) {
  xf <- floor(x)
  df <- x - xf
  if (df < cn & x >= 1) {
    val <- xf - 0.5 + df / (2 * cn)
  }
  if (any(cn <= df & df < (1 - cn), x < 1)) {
    val <- xf
  }

  if (df >= (1 - cn)) {
    val <- xf + 0.5 + (df - 1) / (2 * cn)
  }

  return(val)
}
