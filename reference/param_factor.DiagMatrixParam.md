# Factor of a Diagonal Parameter

Returns \\L = \mathrm{diag}(\sqrt{h(\eta_1)}, \dots,
\sqrt{h(\eta_p)})\\. The Cholesky factor of a diagonal matrix is the
diagonal of its square roots, so this is \\p\\ square roots and no
factorization, against the base class's \\O(p^3)\\.

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` diagonal numeric matrix with a positive
diagonal, satisfying `L %*% t(L) == param_value(s, eta)`, and carrying
no dimnames.

## See also

[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the generic and
[`param_factor.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_factor.matrix_parameter.md)
for what a family without a closed form pays.
