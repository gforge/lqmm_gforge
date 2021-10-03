quad <- function(q, k, type = c("normal", "robust"), rule = 1) {
  if (!(rule %in% 1:4)) {
    warning(paste("Rule ", rule, " not recognised. Rule 1 used (see details in '?lqmm'", sep = ""))
    rule <- 1
  }

  if (rule == 1) {
    odd <- (k %% 2) > 0
    if (type == "normal") {
      QUAD <- gauss.quad.prob(k, type)
      if (odd) QUAD$nodes[floor(k / 2) + 1] <- 0 # ensure middle value is 0
      QUAD$nodes <- permutations(n = k, r = q, v = QUAD$nodes, set = FALSE, repeats.allowed = TRUE)
      QUAD$weights <- apply(permutations(n = k, r = q, v = QUAD$weights, set = FALSE, repeats.allowed = TRUE), 1, prod)
    }

    if (type == "robust") {
      QUAD <- createLaguerre(rule = "product", q = q, k = k)
    }
  } else if (rule == 2) {
    if (k > 25) stop("This rule is for k < 25 only")

    if (type == "normal") {
      QUAD <- createSparseGrid(type = "GQN", dimension = q, k = k)
    }

    if (type == "robust") {
      QUAD <- createLaguerre(rule = "sparse", q = q, k = k)
    }
  } else if (rule == 3) {
    if (k > 25) stop("This rule is for k < 25 only")

    QUAD <- createSparseGrid(type = "KPN", dimension = q, k = k)

    if (type == "robust") {
      warning("Nested rule for integral with Laguerre weights not implemented. Gaussian weights used")
    }
  } else {
    if (k > 25) stop("This rule is for k < 25 only")
    rule <- 2

    if (type == "normal") {
      QUAD <- createSparseGrid(type = "GQN", dimension = q, k = k)
      if (length(QUAD$weights) > k^q) {
        QUAD <- createProductRuleGrid(type = "GQN", dimension = q, k = k)
        rule <- 1
      }
    }

    if (type == "robust") {
      QUAD <- createLaguerre(rule = "sparse", q = q, k = k)
      if (length(QUAD$weights) > k^q) {
        QUAD <- createLaguerre(rule = "product", q = q, k = k)
        rule <- 1
      }
    }
  }

  attr(QUAD, "rule") <- rule
  attr(QUAD, "dimension") <- q
  attr(QUAD, "k") <- k

  return(QUAD)
}
