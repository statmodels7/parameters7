# Third Derivatives of a Scaled Parameter

Closed form: \\\partial^3\_\eta M = h'''(\eta)\\P\\, from
[`linkfunctions7::d3linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkinv.html).
Every order is the same fixed matrix times a link derivative, so this
family costs nothing at any order and is exact at all of them. Under the
default log link it is the matrix again.

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

A list with one matrix keyed as `param_tuple_names(s, 3)`, or an empty
list for a fixed parameter.

## See also

[`param_d4.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.ScaledMatrixParam.md)
for the order above.
