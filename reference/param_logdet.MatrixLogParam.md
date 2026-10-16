# Log-Determinant of a Matrix Logarithm Parameter

Returns the log-determinant in closed form; it is linear in the free
vector. The eigenvalues of \\M = \exp(S)\\ are \\e^{\lambda_i}\\, so

\$\$\log\|M\| = \sum_i \lambda_i = \mathrm{tr}(S) = \sum\_{i=1}^{p}
\eta_i,\$\$

which is the sum of the diagonal free values. The computation is one
[`sum()`](https://rdrr.io/r/base/sum.html) over \\p\\ numbers, without
an eigendecomposition or a determinant, and its cost does not grow in
\\p\\ beyond the sum itself.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number.

## See also

[`param_dlogdet.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.MatrixLogParam.md)
for its gradient, and
[`param_logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.LogCholeskyParam.md),
which is linear for the same kind of reason.
