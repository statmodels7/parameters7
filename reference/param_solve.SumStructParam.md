# Solve and Factor of a Sum of Fixed Matrices

Both come from the assembled matrix: a sum of fixed matrices has no
structure a solve could exploit, unlike the families whose factor is
written out.

## Arguments

- s:

  A
  [`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- b:

  A matrix with `s@dimension` rows, or `NULL`.

- ...:

  Unused.

## Value

A matrix.
