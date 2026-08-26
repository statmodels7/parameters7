# Log-Determinant Hessian of a Diagonal Parameter

Closed form. The log-determinant is a sum of \\\log h(\eta_k)\\ terms,
one free value each, so every mixed component is exactly zero and a pure
one is

\$\$\partial\_{kk} \log\|M\| = m_k\left\\\frac{h''(\eta_k)}{h(\eta_k)} -
\left(\frac{h'(\eta_k)}{h(\eta_k)}\right)^{2}\right\\,\$\$

the second derivative of \\\log h\\ times the number of entries the
value owns.

Under the **log** link this is identically zero, \\\log h(\eta) = \eta\\
being linear, so a check of this quantity against a numerical reference
on a default
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
compares two zeros and would pass whatever was missing. Under a
square-root link at \\\eta = (1, 2)\\ it is \\(-2, -0.5)\\, which is
where the formula has content.

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

A numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
`param_tuple_names(s)` and in that order, the mixed ones exactly zero.

## See also

[`param_dlogdet.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DiagMatrixParam.md)
for the order below, and
[`diag_logdet_higher()`](https://statmodels7.github.io/parameters7/reference/diag_logdet_higher.md)
for orders three and four.
