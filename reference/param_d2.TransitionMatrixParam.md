# Second Derivatives of a Transition Matrix Parameter

Closed form, and mostly zero. A component pairing free values from two
**different** rows is exactly the zero matrix, the rows being
parametrized independently; one pairing two free values of the same row
is that row's
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
second derivative, embedded in that row.

Measured at \\K = 3\\: `alr1.1:alr2.1` is 0 exactly and `alr1.1:alr1.1`
is not. Of the 21 components at \\K = 3\\, only 9 can be non-zero.

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

A list of `choose(s@n_free + 1, 2)` \\K \times K\\ matrices keyed as
`param_tuple_names(s)` and in that order, the cross-row ones exactly
zero.

## See also

[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which assembles it, and
[`param_d1.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.TransitionMatrixParam.md)
and
[`param_d3.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.TransitionMatrixParam.md)
for the neighboring orders.
