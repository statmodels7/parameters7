# Log-Determinant Gradient of a Matrix Logarithm Parameter

Closed form and constant: 1 in each of the \\p\\ diagonal directions and
0 in the \\p(p-1)/2\\ below-diagonal ones, the log-determinant being
\\\sum_i \eta_i\\. The value does not depend on `eta`, which is read
only for its length.

The gradient is 1 where that of
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
is 2: there the diagonal free value is \\\log L\_{ii}\\ and \\\|M\| =
\|L\|^2\\, while here it is a diagonal entry of \\S\\, and \\\log\|M\| =
\mathrm{tr}(S)\\ contains it once.

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

A numeric vector of length `s@n_free`, named by `s@free_names`, holding
1 and 0.

## See also

[`param_logdet.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.MatrixLogParam.md)
for the quantity differentiated, and
[`param_dlogdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.LogCholeskyParam.md),
which returns 2.
