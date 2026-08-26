# Log-Determinant Hessian of a Log-Cholesky Parameter

Closed form, and identically zero: \\\log\|M\| = 2\sum_i \eta_i\\ is
linear in the free vector, so its second derivative vanishes at every
\\\eta\\. The zeros are exact, so a consumer can drop the term entirely
instead of carrying small numbers through a contraction.

This is worth knowing when checking another family: a comparison of a
log-determinant Hessian against a numerical reference can pass on
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
while a term is missing, both sides being zero.
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
are the families where this order has something to get wrong.

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
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

[`param_dlogdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.LogCholeskyParam.md)
for the order below, and
[`param_d3logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.LogCholeskyParam.md),
zero for the same reason.
