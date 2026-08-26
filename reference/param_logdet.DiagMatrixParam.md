# Log-Determinant of a Diagonal Parameter

Returns \\\log\|M\| = \sum\_{i=1}^{p} \log h(\eta_i)\\, the sum of the
logarithms of the diagonal entries, which for a diagonal matrix is the
log-determinant exactly. For a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
the one entry is counted \\p\\ times, giving \\p \log h(\eta_1)\\.

Under the **log** link \\\log h(\eta) = \eta\\, so the result is
`sum(eta)` and the quantity is linear in the free vector; under any
other link it is not.

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

A single number.

## See also

[`param_dlogdet.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DiagMatrixParam.md)
for its gradient, and
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the derivatives of \\\log h\\.
