# Derivatives of a Block-Diagonal Parameter

Each block's own derivatives, placed in the rows and columns that block
occupies. A component whose indices span two blocks is exactly zero, the
free values of one block not entering another.

## Arguments

- s:

  A
  [`BlockDiagParam`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Unused.

## Value

A named list of symmetric matrices.
