# Log-Determinant of a Block-Diagonal Parameter

Returns the sum of the blocks' log-(pseudo-)determinants, because the
determinant of a block-diagonal matrix is the product of the blocks'
determinants. Each block computes its own by whatever route is available
to it, so a closed form in a block stays a closed form in the composite,
and the method does not take the determinant of the assembled matrix.

Where a block is rank deficient its term is the log pseudo-determinant,
the sum over the non-zero eigenvalues, and the composite's is then a
pseudo one too.

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

A single number.

## See also

[`param_dlogdet.BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.BlockDiagParam.md)
for its derivatives.
