# Default Factor

The method every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
of its own. It evaluates
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and returns the lower triangular Cholesky factor through
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md),
which decides positive definiteness from the eigenvalues before
attempting the factorization. Exact, at \\O(p^3)\\.

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object, of full rank.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` lower triangular numeric matrix with a
positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)`.

## Details

A family that declares full rank and is then not positive definite at
this \\\eta\\ throws, naming the family and saying that the verdict is
spectral. That distinction matters to whoever reads the message: a
caught [`chol()`](https://rdrr.io/r/base/chol.html) error would be a
statement about the arithmetic, and could differ between platforms on a
matrix with an exactly zero eigenvalue, while a test on the eigenvalues
is a statement about the matrix. The error means the family's own
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
has left the cone it claims to parametrize.

## See also

[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the generic, which rejects a rank-deficient family before this runs,
and
[`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
for the definiteness test.
