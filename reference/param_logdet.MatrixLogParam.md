# Log-Determinant of a Matrix Logarithm Parameter

Closed form and linear. The eigenvalues of \\M = \exp(S)\\ are
\\e^{\lambda_i}\\, so

\$\$\log\|M\| = \sum_i \lambda_i = \mathrm{tr}(S) = \sum\_{i=1}^{p}
\eta_i,\$\$

the sum of the diagonal free values. One
[`sum()`](https://rdrr.io/r/base/sum.html) over \\p\\ numbers: no
eigendecomposition and no determinant, with nothing growing in \\p\\
beyond the sum itself. It agrees with the eigenvalues of the assembled
matrix to the printed digit.

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
