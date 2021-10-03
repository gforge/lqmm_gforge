## Use `SparseGrid' version 0.8.1 (Jelmer Ypma)

createLaguerre <- function(rule, q, k) {
  laguerre <- function(k) {
    odd <- (k %% 2) > 0
    if (odd) {
      val <- gauss.quad(n = (k - 1) / 2, kind = "laguerre")
      val$nodes <- c(-rev(val$nodes), 0, val$nodes)
      val$weights <- c(rev(val$weights), 0, val$weights)
    } else {
      val <- gauss.quad(n = k / 2, kind = "laguerre")
      val$nodes <- c(-rev(val$nodes), val$nodes)
      val$weights <- c(rev(val$weights), val$weights)
    }
    return(val)
  }

  if (rule == "product") {
    QUAD <- laguerre(k)
    QUAD$nodes <- permutations(n = k, r = q, v = QUAD$nodes, set = FALSE, repeats.allowed = TRUE)
    QUAD$weights <- apply(permutations(n = k, r = q, v = QUAD$weights, set = FALSE, repeats.allowed = TRUE), 1, prod)
  }


  if (rule == "sparse") {
    QUAD <- suppressWarnings(createSparseGrid(laguerre, dimension = q, k = k, sym = TRUE))
    QUAD$weights <- QUAD$weights * 2
  }

  return(QUAD)
}
