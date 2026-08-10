# Free Vector of a Scales-Times-Correlation Parameter

Reads the standard deviations off the diagonal, divides them out, and
hands the resulting correlation matrix to the correlation block.

## Arguments

- s:

  A
  [`DrProdParam`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- m:

  A symmetric positive definite matrix of side `s@dimension`.

- ...:

  Unused.

## Value

A numeric vector of length `s@n_free`.
