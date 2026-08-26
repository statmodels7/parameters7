# Value of a Transition Matrix Parameter

Returns the transition matrix: each row is the softmax of that row's own
\\K-1\\ free values, through
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md),
so each row is positive and sums to exactly 1. Nothing is tested and
nothing is renormalized; the row sums are a property of the map. Rows
are the distributions, so \\P\_{ij}\\ is the probability of moving from
state \\i\\ to state \\j\\.

## Arguments

- s:

  A
  [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A \\K \times K\\ numeric matrix with positive entries and rows summing
to 1, with dimnames `s1`, `s2`, ... on both margins.

## See also

[`param_free.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_free.TransitionMatrixParam.md)
for the inverse, and
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
for the per-row arithmetic.
