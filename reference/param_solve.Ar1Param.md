# Solve of an AR(1) Parameter

Exact and tridiagonal. The precision of an AR(1) pattern is
\\(1-\rho^2)^{-1}\\ times the matrix with \\1\\ at the two corners of
the diagonal, \\1+\rho^2\\ elsewhere on it and \\-\rho\\ on the first
off-diagonals; no factorisation is performed.

## Arguments

- s:

  An
  [`Ar1Param`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values.

- b:

  A numeric matrix with `s@dimension` rows.

- ...:

  Unused.

## Value

A numeric matrix.
