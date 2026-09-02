# Derivative Arrays of an Inverse Parameter

The four orders of \\\partial N\\ for \\N = S^{-1}\\, each the ordered
set-partition sum
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
writes out. Every factor comes from the inner family's own arrays, so a
closed form there stays closed here.

## Arguments

- s:

  An
  [`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A named list of `s@dimension` square matrices, keyed as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
says.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the formula and
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
for the contract.
