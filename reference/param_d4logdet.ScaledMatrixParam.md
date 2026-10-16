# Fourth Log-Determinant Derivatives of a Scaled Parameter

The fourth derivative has closed form: \\r\\ times the fourth derivative
of \\\log h\\, which
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
gives as \\u_4 - 4u_1u_3 - 3u_2^2 + 12u_1^2u_2 - 6u_1^4\\ with \\u_m =
h^{(m)}/h\\. The constant \\\log\|P\|\_+\\ contributes nothing.

The quantity is exactly zero under the default log link, because it is
linear in the free value. This is the order at which a numerical route
is least accurate, so a family with a curved link gains most from the
closed form here.

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

A numeric vector with one entry keyed as `param_tuple_names(s, 4)`, or
empty for a fixed parameter.

## See also

[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the expression, and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
for the numerical route.
