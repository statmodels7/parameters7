# Derivatives of a Scales-Times-Correlation Parameter

Each component is the scale factor times the correlation's own
component, elementwise, the two groups of free values being disjoint.

## Arguments

- s:

  A
  [`DrProdParam`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
