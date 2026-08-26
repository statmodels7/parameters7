# First Derivative of a Scaled Parameter

Closed form: \\\partial\_\eta M = h'(\eta)\\P\\, the fixed matrix times
the link's own first derivative. Under the default log link \\h' = h\\,
so the derivative **is** the matrix, which the example on
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
checks.

A fixed parameter, built with `link = NULL`, has no free value, so the
list is empty.

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

A list with one symmetric matrix named by `s@free_names`, or an empty
list for a fixed parameter.

## See also

[`param_d2.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.ScaledMatrixParam.md)
for the order above, and
[`scaled_scale()`](https://statmodels7.github.io/parameters7/reference/scaled_scale.md)
for the link derivatives.
