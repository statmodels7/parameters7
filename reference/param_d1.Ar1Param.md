# Derivatives of an AR(1) Parameter

Closed form at every order: a component with \\a\\ scale indices and
\\b\\ correlation indices is the \\a\\-th derivative of the scale times
the \\b\\-th derivative of the pattern, and the pattern's derivatives
are powers composed with the link.

## Arguments

- s:

  An
  [`Ar1Param`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
