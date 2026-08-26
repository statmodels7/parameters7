# Log-Determinant of a Log-Cholesky Parameter

Closed form, and linear in the free vector. Since \\\|M\| = \|L\|^2\\
and \\L\\ is triangular,

\$\$\log\|M\| = 2\sum\_{i=1}^{p} \log L\_{ii} = 2\sum\_{i=1}^{p}
\eta_i,\$\$

twice the sum of the free values on the diagonal. It is one
[`sum()`](https://rdrr.io/r/base/sum.html) over \\p\\ numbers: no
factorization, no determinant, no eigendecomposition, and nothing that
grows with \\p\\ beyond the sum itself, against the base class's
\\O(p^3)\\.

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number, always finite: `eta` is finite by the generic's check
and the sum of \\p\\ finite numbers cannot overflow at any usable \\p\\.

## See also

[`param_dlogdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.LogCholeskyParam.md)
for its gradient, and
[`param_logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_logdet.matrix_parameter.md)
for the route a family without a closed form takes.
