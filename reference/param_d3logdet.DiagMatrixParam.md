# Third Log-Determinant Derivatives of a Diagonal Parameter

Closed form. The log-determinant is a sum of \\\log h(\eta_k)\\ terms,
so each pure component is the third derivative of \\\log h\\ at that
free value, times the number of entries the value owns, and every mixed
component is exactly zero.

Under the log link the whole vector is zero, \\\log h(\eta) = \eta\\
being linear. A link with curvature is where this order has content.

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

A numeric vector of `choose(s@n_free + 2, 3)` entries, keyed as
`param_tuple_names(s, 3)`.

## See also

[`diag_logdet_higher()`](https://statmodels7.github.io/parameters7/reference/diag_logdet_higher.md),
which assembles it, and
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the derivatives of \\\log h\\.
