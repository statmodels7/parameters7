# Derivative Arrays of an Inverse Autoregressive Parameter

The four orders of \\\partial\Omega\\ for \\\Omega\\ the precision of an
autoregression of order \\q\\, written out from the prediction
factorization instead of assembled by the ordered-block-partition sum
that
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
uses. The two routes agree to rounding.

## Arguments

- s:

  An
  [`AutoregressiveInvParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveInvParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A named list of `s@dimension` square matrices, keyed as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
gives them.

## See also

[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
for the formulas.
