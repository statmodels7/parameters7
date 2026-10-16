# Log-Determinant Hessian of a Log-Cholesky Parameter

Computes the second derivatives in closed form, and they are identically
zero: \\\log\|M\| = 2\sum_i \eta_i\\ is linear in the free vector, so
its second derivative vanishes at every \\\eta\\. The zeros are exact,
so a consumer can drop the term entirely instead of carrying small
numbers through a contraction.

Therefore, when a general check compares the log-determinant Hessian
with a numerical reference on
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
both sides are zero.
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
and
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
are families where this order is not zero.

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
