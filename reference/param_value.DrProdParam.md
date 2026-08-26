# Value of a Scales-Times-Correlation Parameter

Forms \\D R D\\ as `outer(d, d) * R`, the elementwise product being what
two diagonal multiplications amount to. The correlation block is asked
for its own value at its own stretch of the free vector, so the
composite is exactly what that family returns, rescaled.

The value is labeled `v1`, `v2`, ..., `vp` on both margins, the
convention
[`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
states and every family in the package follows.

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A symmetric positive definite `s@dimension` by `s@dimension` numeric
matrix, labeled `v1`, `v2`, ..., `vp` on both margins.

## See also

[`param_free.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_free.DrProdParam.md)
for the inverse, and
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
for the parametrization.
