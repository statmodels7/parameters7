# Default Factor

Fallback: the lower Cholesky factor of
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md),
refused when the matrix is not positive definite spectrally.

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Unused.

## Value

A lower triangular numeric matrix.
