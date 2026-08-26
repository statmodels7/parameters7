# Log-Determinant Gradient of a Block Replication

Closed form: \\m\\ times the inner parameter's gradient. A shared free
value moves every block, so the gradient is multiplied rather than
copied: at \\m = 3\\ over a 2 x 2 log-Cholesky covariance it is
`c(6, 6, 0)` where the inner parameter answers `c(2, 2, 0)`.

## Arguments

- s:

  A
  [`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## See also

[`param_logdet.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.KronIdentityParam.md)
for the quantity differentiated.
