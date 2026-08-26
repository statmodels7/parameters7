# Fourth Derivatives of a Scaled Parameter

Closed form: \\\partial^4\_\eta M = h''''(\eta)\\P\\, from
[`linkfunctions7::d4linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkinv.html).
Exact, where a family without a closed form would get a product stencil
good to about five digits at this order.

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

A list with one matrix keyed as `param_tuple_names(s, 4)`, or an empty
list for a fixed parameter.

## See also

[`param_d3.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.ScaledMatrixParam.md)
for the order below, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the alternative.
