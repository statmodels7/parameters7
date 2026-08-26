# Second Derivatives of a Matrix Logarithm Parameter

The second Frechet derivative: chains of two rotated directions
contracted against three-point divided differences \\e\[\lambda_i,
\lambda_j, \lambda_k\]\\, summed over both orderings of the two
directions.

Summing over the orderings, never over the distinct ones, is what brings
a repeated index out right: the component in one free value twice gets
its factor of 2 from the sum instead of from a correction applied
afterwards.

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

A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
`param_tuple_names(s)` and in that order.

## See also

[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
which does the contraction, and
[`param_d1.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.MatrixLogParam.md)
and
[`param_d3.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.MatrixLogParam.md)
for the neighboring orders.
