# Value and Inverse Map of an Inverse Parameter

[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns the inner family's inverse, read through
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
so that a family with a closed inverse never reaches a factorization.
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
inverts the matrix given and hands it to the inner family's own inverse
map, so the free vector is the inner one and the round trip closes on
it.

## Arguments

- s:

  An
  [`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- m:

  A symmetric positive definite matrix of side `s@dimension`.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
a `s@dimension` square matrix,
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
a numeric vector of length `s@n_free`.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the family.
