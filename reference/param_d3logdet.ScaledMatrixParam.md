# Third Log-Determinant Derivatives of a Scaled Parameter

Closed form. The log-(pseudo-)determinant is \\r \log h(\eta) +
\log\|P\|\_+\\ with \\r\\ the rank, so every derivative is \\r\\ times
the matching derivative of \\\log h\\, and the constant contributes
nothing. With one free value there is one component.

Under the default log link \\\log h(\eta) = \eta\\, so the answer is
exactly zero and this family cannot exercise the order; a link with
curvature can.

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- eta:

  A numeric vector with one free value, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector with one entry keyed as `param_tuple_names(s, 3)`, or
empty for a fixed parameter.

## See also

[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md),
which supplies the derivative of \\\log h\\, and
[`param_d4logdet.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.ScaledMatrixParam.md)
for the order above.
