# Second Derivative of a Scaled Parameter

Closed form: \\\partial^2\_\eta M = h''(\eta)\\P\\. With one free value
there is one component, and under the default log link it is the matrix
again, \\h'' = h\\. Empty for a fixed parameter.

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

A list with one symmetric matrix keyed as `param_tuple_names(s)`, or an
empty list for a fixed parameter.

## See also

[`param_d1.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.ScaledMatrixParam.md)
and
[`param_d3.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.ScaledMatrixParam.md)
for the neighboring orders.
