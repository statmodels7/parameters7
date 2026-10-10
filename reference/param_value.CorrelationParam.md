# Value of a Correlation Parameter

Returns \\R = LL^\top\\, with \\L\\ assembled from the angles: row \\i\\
is a unit vector in spherical coordinates, so the diagonal of \\R\\ is 1
and its off-diagonal entries are correlations. The diagonal is then set
to exactly 1, which removes the rounding of the telescoping sum; the
positive definiteness follows from the construction and is not tested.
The cost is one \\O(p^3)\\ product.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A `s@dimension` by `s@dimension` correlation matrix: symmetric, positive
definite, unit diagonal, dimnames `v1`, `v2`, ...

## See also

[`param_free.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_free.CorrelationParam.md)
for the inverse, and
[`param_factor.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_factor.CorrelationParam.md)
for \\L\\ itself.
