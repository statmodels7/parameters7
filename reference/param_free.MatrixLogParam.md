# Free Vector of a Matrix Logarithm Parameter

Returns the matrix logarithm of `m`, read off its eigendecomposition as
\\Q \log(\Lambda) Q^\top\\ and then read out of the lower triangle.
Exact and a true inverse of
[`param_value.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_value.MatrixLogParam.md),
the matrix logarithm of a symmetric positive definite matrix being
unique among symmetric matrices. The round trip passes through two
eigendecompositions, one in each direction.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- m:

  A symmetric positive definite `s@dimension` by `s@dimension` numeric
  matrix, already checked for shape and symmetry by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

A matrix that is not positive definite is rejected: the logarithm of a
non-positive eigenvalue is not a real number, so there is no free vector
to return. The test is on the eigenvalues of `m`.

## See also

[`param_value.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_value.MatrixLogParam.md),
the map this inverts.
