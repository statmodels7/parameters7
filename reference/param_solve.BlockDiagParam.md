# Solve and Factor of a Block-Diagonal Parameter

Both are blockwise: the inverse of a block-diagonal matrix is the block
diagonal of the inverses, and the same holds of a triangular factor.

## Arguments

- s:

  A
  [`BlockDiagParam`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- b:

  A matrix with `s@dimension` rows, or `NULL`.

- ...:

  Unused.

## Value

A matrix.
