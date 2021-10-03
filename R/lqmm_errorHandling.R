errorHandling <- function(code, type, maxit, tol, fn) {
  txt <- switch(type,
    low = "Lower loop",
    upp = "Upper loop"
  )

  if (code == -1) {
    warning(paste(txt, " did not converge in: ", fn, ". Try increasing max number of iterations ", "(", maxit, ") or tolerance (", tol,
      ")\n",
      sep = ""
    ))
  }
  if (code == -2) warning(paste(txt, " did not start in: ", fn, ". Check max number of iterations ", "(", maxit, ")\n", sep = ""))
}
