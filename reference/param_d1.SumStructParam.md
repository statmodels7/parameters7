# Derivatives of a Sum of Fixed Matrices

The value being linear in the weights, a component is zero unless every
index names the same free value, and is then that weight's derivative
times its component.

## Arguments

- s:

  A
  [`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
