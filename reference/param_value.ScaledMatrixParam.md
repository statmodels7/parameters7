# Matrix of a Scaled Parameter

Returns \\M = h(\eta)\\P\\, the stored fixed matrix times the scale,
which is \\P\\ itself for a fixed parameter. The computation is one
elementwise multiplication and does not decompose the matrix. The result
inherits \\P\\'s rank exactly, because a positive multiple leaves a null
space alone, and the class can therefore record that rank at
construction.

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A symmetric positive semidefinite `s@dimension` by `s@dimension` numeric
matrix with dimnames `v1`, `v2`, ..., of rank `s@rank`. The matrix is
positive definite only where `s@rank == s@dimension`.

## See also

[`param_free.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_free.ScaledMatrixParam.md)
for the inverse, and
[`scaled_scale()`](https://statmodels7.github.io/parameters7/reference/scaled_scale.md)
for the scale.
