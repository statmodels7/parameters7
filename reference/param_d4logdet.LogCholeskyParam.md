# Fourth Log-Determinant Derivatives of a Log-Cholesky Parameter

Computes the fourth derivatives in closed form, and they are identically
zero, for the reason that they are zero at second and third order:
\\\log\|M\| = 2\sum_i \eta_i\\ is linear in the free vector. The zeros
are exact.

A numerical fallback is least accurate at this order, so a closed form
gains the most here for a family whose log-determinant is not linear.
See
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the fallback.

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

A numeric vector of `choose(s@n_free + 3, 4)` zeros, keyed as
`param_tuple_names(s, 4)`.

## See also

[`param_d3logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.LogCholeskyParam.md)
for the order below, and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the numerical route.
