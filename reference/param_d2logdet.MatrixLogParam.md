# Log-Determinant Hessian of a Matrix Logarithm Parameter

Closed form, and identically zero: \\\log\|M\| = \mathrm{tr}(S)\\ is
linear in the free vector, so its second derivative vanishes at every
\\\eta\\. The zeros are exact.

As with
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
that means this family cannot exercise the order: a check of a
log-determinant Hessian against a numerical reference passes here
whatever is missing.
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
are where it has content.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic. Its values do not enter the result.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of `choose(s@n_free + 1, 2)` zeros, keyed as
`param_tuple_names(s)`.

## See also

[`param_dlogdet.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.MatrixLogParam.md)
for the order below, and
[`param_d3logdet.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.MatrixLogParam.md),
zero for the same reason.
