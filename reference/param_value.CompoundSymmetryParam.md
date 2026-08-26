# Value of a Compound Symmetry Parameter

Returns \\\sigma^2\\(1-\rho)I + \rho J\\\\: the common variance on the
diagonal and the common covariance \\\sigma^2\rho\\ everywhere else.
Positive definiteness follows from the correlation's link, which is
bounded below at \\-1/(p-1)\\, so nothing is tested here.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` symmetric positive definite numeric
matrix with a constant diagonal and constant off-diagonal entries, with
dimnames `v1`, `v2`, ...

## See also

[`param_free.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_free.CompoundSymmetryParam.md)
for the inverse map, and
[`param_solve.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_solve.CompoundSymmetryParam.md)
for the closed matrix inverse.
