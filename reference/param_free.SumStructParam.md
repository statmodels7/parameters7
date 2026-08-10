# Free Vector of a Sum of Fixed Matrices

Recovers the weights by least squares on the components' entries and
rejects a matrix the combination cannot reproduce. Non-positive weights
are rejected too, the family carrying them through a positive link.

## Arguments

- s:

  A
  [`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- m:

  A symmetric matrix of side `s@dimension`.

- ...:

  Unused.

## Value

A numeric vector of length `s@n_free`.
