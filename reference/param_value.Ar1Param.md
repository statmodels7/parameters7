# Value of an AR(1) Parameter

Returns \\M\_{ij} = \sigma^2 \rho^{\|i-j\|}\\: the common variance on
the diagonal and a covariance falling geometrically with the lag.
Positive definiteness holds at every free vector, the rhobit link
keeping \\\|\rho\| \< 1\\, so nothing is tested here.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` symmetric positive definite Toeplitz
matrix with dimnames `v1`, `v2`, ...

## See also

[`param_free.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_free.Ar1Param.md)
for the inverse map,
[`param_solve.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_solve.Ar1Param.md)
for the tridiagonal matrix inverse, and
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
for the pattern.
