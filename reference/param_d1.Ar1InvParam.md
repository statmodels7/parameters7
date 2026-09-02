# Derivative Arrays of an Inverse AR(1) Parameter

The four orders of \\\partial\Omega\\ for \\\Omega\\ the precision of an
AR(1), written out from the product structure rather than assembled by
the ordered-block-partition sum
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
uses. The two routes agree to machine precision.

## Arguments

- s:

  An
  [`Ar1InvParam()`](https://statmodels7.github.io/parameters7/reference/Ar1InvParam.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A named list of `s@dimension` square matrices, keyed as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
says.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
for the formulas.
