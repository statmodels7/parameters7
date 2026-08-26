# Free Vector of a Transition Matrix Parameter

Returns the additive log-ratio of each row, \\\log(P\_{ij}/P\_{iK})\\,
exact and a true inverse of
[`param_value.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.TransitionMatrixParam.md):
the round trip closes to \\2 \times 10^{-16}\\.

## Arguments

- s:

  A
  [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  object.

- m:

  A \\K \times K\\ row-stochastic numeric matrix: strictly positive with
  every row summing to 1.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

A row is rejected when it has a non-positive entry, being then outside
the **open** simplex where \\\log 0\\ is not finite, or when it does not
sum to 1, in which case it is not a probability distribution and is
**not** renormalized. The message names the row.

The first rejection is the one a fit meets:
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
saturates a large free value to an exact 0 in the reference state, so a
matrix produced at the boundary cannot be inverted back.

## See also

[`param_value.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.TransitionMatrixParam.md),
the map this inverts, and
[`param_free.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_free.SimplexParam.md),
which does one row's worth.
