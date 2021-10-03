switch_check <- function(l0, l1, tol_ll, t0, t1, tol_theta, rule = "1") {
  deltal <- abs(l1 / l0 - 1)
  deltat <- max(abs(t1 / t0 - 1))

  switch(rule,
    "1" = deltal < tol_ll,
    "2" = deltal < tol_ll && deltat < tol_theta
  )
}
