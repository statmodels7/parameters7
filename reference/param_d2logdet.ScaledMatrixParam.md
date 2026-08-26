# Log-Determinant Hessian of a Scaled Parameter

Closed form: \\r\\h''/h - (h'/h)^2\\\\, the rank times the second
derivative of \\\log h\\. Under the default log link it is exactly zero,
the log pseudo-determinant being \\r\eta + \log\|P\|\_+\\ and so linear
in the free value. Under a link with curvature it is not.

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector with one entry keyed as `param_tuple_names(s)`, or
empty for a fixed parameter.

## See also

[`param_dlogdet.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.ScaledMatrixParam.md)
for the order below, and
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for orders three and four.
