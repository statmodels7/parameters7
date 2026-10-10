# Fourth Derivatives of a Matrix Logarithm Parameter

Chains of four rotated directions contracted against five-point divided
differences, summed over the orderings of the four directions, up to
twenty-four per component. Exact, and much more expensive than the
fourth derivatives of
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
on a matrix of the same size.

The cost is the price of a linear log-determinant and of an inverse
equal to the map at \\-\eta\\. When fourth derivatives are taken
repeatedly and those two properties are not needed,
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
parametrizes the same cone at a much lower cost.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
`param_tuple_names(s, 4)` and in that order.

## See also

[`mlog_higher()`](https://statmodels7.github.io/parameters7/reference/mlog_higher.md),
which assembles it,
[`param_d3.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.MatrixLogParam.md)
for the order below, and
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for the cheaper chart.
