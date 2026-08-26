# Log-Determinant Gradient of a Log-Cholesky Parameter

Closed form and constant: 2 in each of the \\p\\ diagonal directions and
0 in the \\p(p-1)/2\\ below-diagonal ones, since \\\log\|M\|\\ is
\\2\sum_i \eta_i\\ and does not involve the rest of \\\eta\\ at all. The
value does not depend on `eta`, which is read only for its length.

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

A numeric vector of length `s@n_free`, named by `s@free_names`, holding
2 and 0.

## See also

[`param_logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.LogCholeskyParam.md)
for the quantity differentiated, and
[`param_d2logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.LogCholeskyParam.md),
which is zero for the same reason.
