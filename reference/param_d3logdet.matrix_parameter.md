# Default Higher Log-Determinant Derivatives

Fallback: one central stencil on
[`param_d2logdet`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
which is the exact trace identity given the matrix derivatives – a
single layer on an analytic quantity, per the toolkit's rule.

## Arguments

- s:

  A
  [`matrix_parameter`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named numeric vector.
