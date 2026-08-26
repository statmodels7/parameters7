# Third Log-Determinant Derivatives of a Matrix Logarithm Parameter

Closed form, and identically zero, the trace being linear in the free
vector. The zeros are exact, so a consumer can drop the term instead of
carrying small numbers through a contraction.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic. Its values do not enter the result.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of `choose(s@n_free + 2, 3)` zeros, keyed as
`param_tuple_names(s, 3)`.

## See also

[`param_d2logdet.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.MatrixLogParam.md)
for the order below, and
[`param_d4logdet.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.MatrixLogParam.md)
for the one above.
