# Solve and Factor of a Scales-Times-Correlation Parameter

\\\Sigma^{-1} = D^{-1} R^{-1} D^{-1}\\ and \\\Sigma = (DL)(DL)'\\ for
\\L\\ the correlation's factor, so both come from the correlation block
with a scaling on either side.

## Arguments

- s:

  A
  [`DrProdParam`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- b:

  A matrix with `s@dimension` rows, or `NULL`.

- ...:

  Unused.

## Value

A matrix.
