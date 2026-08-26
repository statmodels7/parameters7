# Free Vector of a Log-Cholesky Parameter

Returns the free vector behind a matrix: the Cholesky factor of `m`,
read off at the positions
[`chol_positions()`](https://statmodels7.github.io/parameters7/reference/chol_positions.md)
records, with the diagonal logged. Exact, and a true inverse of
[`param_value.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_value.LogCholeskyParam.md)
because the triangular factor with a positive diagonal is unique.
Measured, the round trip closes to \\7 \times 10^{-17}\\.

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- m:

  A symmetric positive definite `s@dimension` by `s@dimension` numeric
  matrix, already checked for shape and symmetry by the generic. A
  matrix that is not positive definite throws.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

`m` is rejected when it is not positive definite, with the verdict taken
from the eigenvalues through
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
and the message saying so. That distinction matters: a caught
[`base::chol()`](https://rdrr.io/r/base/chol.html) error would be a
statement about the arithmetic and can differ between platforms on a
matrix with an exactly zero eigenvalue, while a test on the spectrum is
a statement about the matrix.

The logarithm is applied by subsetting on `on_diagonal`. Through
[`ifelse()`](https://rdrr.io/r/base/ifelse.html) it would be evaluated
over the whole vector, including the below-diagonal entries that are
free to be negative, which produces `NaN`s and a warning about values
the function then discards.

## See also

[`param_value.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_value.LogCholeskyParam.md),
the map this inverts, and
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
for the definiteness test.
