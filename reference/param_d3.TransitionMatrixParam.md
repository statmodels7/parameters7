# Third Derivatives of a Transition Matrix Parameter

Closed form, from
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)'s
cumulant recursion at third order, embedded row by row. A component
whose three free values do not all belong to one row is exactly the zero
matrix, so the great majority of the list is zero and is skipped instead
of computed.

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

A list of `choose(s@n_free + 2, 3)` \\K \times K\\ matrices keyed as
`param_tuple_names(s, 3)` and in that order.

## See also

[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which assembles it, and
[`param_d4.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.TransitionMatrixParam.md)
for the order above.
