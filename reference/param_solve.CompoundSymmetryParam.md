# Solve of a Compound Symmetry Parameter

Exact, by Sherman-Morrison: the inverse of \\\sigma^2\\(1-\rho)I + \rho
J\\\\ is \\\\\sigma^2(1-\rho)\\^{-1}\[I - \rho J/\\1 + (p-1)\rho\\\]\\,
compound symmetric again, so no factorization is performed.

## Arguments

- s:

  A
  [`CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object.

- eta:

  A numeric vector of two free values.

- b:

  A numeric matrix with `s@dimension` rows.

- ...:

  Unused.

## Value

A numeric matrix.
