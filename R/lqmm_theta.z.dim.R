theta.z.dim <- function(type, n) {
  switch(type,
    "pdIdent" = 1,
    "pdDiag" = n,
    "pdCompSymm" = if (n == 1) 1 else 2,
    "pdSymm" = n * (n + 1) / 2
  )
}
