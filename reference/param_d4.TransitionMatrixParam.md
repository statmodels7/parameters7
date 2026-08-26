# Fourth Derivatives of a Transition Matrix Parameter

Closed form, from
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)'s
cumulant recursion at fourth order, embedded row by row, and the top of
the contract. A component whose four free values do not all belong to
one row is exactly the zero matrix.

Exactness matters most here: a product stencil at fourth order keeps
about five digits, and the list is large, so a numerical route would be
both slow and poor. The row-wise structure is what keeps it affordable
at all, the cross-row components never being evaluated.

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

A list of `choose(s@n_free + 3, 4)` \\K \times K\\ matrices keyed as
`param_tuple_names(s, 4)` and in that order. At \\K = 3\\ that is 126
components, of which 15 can be non-zero.

## See also

[`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md),
which assembles it,
[`param_d3.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.TransitionMatrixParam.md)
for the order below, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the alternative.
