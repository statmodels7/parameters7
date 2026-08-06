# Derivatives of an Autoregressive Parameter

Closed form at every order. The map from the partial autocorrelations to
the matrix is polynomial, so the derivative arrays propagated through
the Levinson-Durbin recursion give each derivative exactly; nothing is
differenced.

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
