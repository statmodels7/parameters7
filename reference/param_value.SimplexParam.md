# Value of a Simplex Parameter

Returns the probability vector \\\pi\\, the softmax of the free vector
with the reference category's implicit 0 appended, named `p1` ... `pK`.
The entries are non-negative and sum to 1 up to rounding at every free
vector, so nothing is tested and nothing is renormalized. Computed
through
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)'s
log-sum-exp shift, so a large free value saturates instead of
overflowing.

## Arguments

- s:

  A
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  object, whose `param_params$n_cat` is read.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length \\K\\ named `p1` ... `pK`, with non-negative
entries summing to 1.

## See also

[`param_free.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_free.SimplexParam.md)
for the inverse, and
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
for the arithmetic.
