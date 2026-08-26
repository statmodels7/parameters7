# Fourth Derivatives of a Matrix Logarithm Parameter

Chains of four rotated directions contracted against five-point divided
differences, summed over the twenty-four orderings. Exact, and the most
expensive quantity the package computes: measured at \\p = 4\\ one call
costs **3.4 s**, against 0.003 s for
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)'s
fourth derivatives on the same matrix size.

That factor of a thousand is the trade this chart makes, and it is worth
knowing before putting a
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
inside a loop that takes fourth derivatives. It buys a linear
log-determinant and an inverse that is a sign flip; if those are not
what the model needs,
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
parametrizes the same cone.

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
