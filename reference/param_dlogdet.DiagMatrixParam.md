# Log-Determinant Gradient of a Diagonal Parameter

Closed form,

\$\$\partial_k \log\|M\| = m_k \frac{h'(\eta_k)}{h(\eta_k)},\$\$

with \\m_k\\ the number of diagonal entries the \\k\\-th free value
owns: 1 for a
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
and \\p\\ for the single value of a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
which enters the log-determinant once per entry.

Under the log link \\h'/h = 1\\, so the answer is a vector of ones for a
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
and the scalar \\p\\ for a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).

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

A numeric vector of length `s@n_free`, named by `s@free_names`.

## See also

[`diag_multiplicity()`](https://statmodels7.github.io/parameters7/reference/diag_multiplicity.md)
for \\m_k\\, and
[`param_d2logdet.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.DiagMatrixParam.md)
for the order above.
