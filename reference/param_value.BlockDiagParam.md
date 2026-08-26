# Value of a Block-Diagonal Parameter

Evaluates each block at its own stretch of the free vector and writes it
into the rows and columns that block occupies, leaving the off-diagonal
blocks at the zeros the matrix was created with. Each block is exactly
what its own family returns, to the bit.

The value is labeled `v1`, `v2`, ..., `vp` on both margins, the
convention
[`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
states and every family in the package follows.

## Arguments

- s:

  A
  [`BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A symmetric `s@dimension` by `s@dimension` numeric matrix, block
diagonal and labeled `v1`, `v2`, ..., `vp` on both margins.

## See also

[`param_free.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_free.BlockDiagParam.md)
for the inverse, and
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
for the composition.
