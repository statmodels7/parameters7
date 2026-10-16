# Solve of a Matrix Logarithm Parameter

Returns \\M^{-1} b\\ exactly, where \\M^{-1} = \exp(-S)\\ is computed
through the same eigendecomposition as the value.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values.

- b:

  A numeric matrix with `s@dimension` rows.

- ...:

  Unused.

## Value

A numeric matrix.
