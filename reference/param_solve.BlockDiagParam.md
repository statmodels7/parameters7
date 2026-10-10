# Solve and Factor of a Block-Diagonal Parameter

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for a
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
parameter, both of them blockwise: the inverse of a block-diagonal
matrix is the block diagonal of the inverses, and the same holds of a
lower triangular factor. Each block uses its own route, so the
tridiagonal inverse of
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and the assembled factor of
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
are used in the composite.

## Arguments

- s:

  A
  [`BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- b:

  A numeric matrix with `s@dimension` rows, defaulted to the identity by
  the generic.
  [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
  takes no `b`.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
returns a numeric matrix with `s@dimension` rows and as many columns as
`b`;
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
a lower triangular `s@dimension` by `s@dimension` matrix.

## Details

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
takes the rows of `b` that belong to each block and hands them to that
block, so the whole matrix is never inverted.
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
assembles the blocks' factors on the diagonal; the result is lower
triangular with \\M = L L^\top\\, as
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
requires.

Both are rejected by the generic where any block is rank deficient, the
composite having no inverse and no Cholesky factor then.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the two generics.
