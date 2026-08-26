# Fourth Log-Determinant Derivatives of a Diagonal Parameter

Closed form, by the same separability as at third order: each pure
component is the fourth derivative of \\\log h\\ at that free value,
times the number of entries it owns, and every mixed component is
exactly zero. The fourth derivative of \\\log h\\ is \\u_4 - 4u_1u_3 -
3u_2^2 + 12u_1^2u_2 - 6u_1^4\\ with \\u_m = h^{(m)}/h\\, written out in
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md).

Under the log link the whole vector is zero. This is the order where a
numerical route is least usable, so the closed form is worth most here.

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

A numeric vector of `choose(s@n_free + 3, 4)` entries, keyed as
`param_tuple_names(s, 4)`.

## See also

[`diag_logdet_higher()`](https://statmodels7.github.io/parameters7/reference/diag_logdet_higher.md),
which assembles it,
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the formula, and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the numerical route.
