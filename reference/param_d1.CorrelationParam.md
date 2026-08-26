# First Derivatives of a Correlation Parameter

Closed form, by the Leibniz rule on \\R = LL^\top\\,

\$\$\partial_k R = L_k L^\top + L L_k^\top,\$\$

with the factor's derivative \\L_k\\ built from the trigonometric tables
of the angles: differentiating an entry of \\L\\ replaces one sine or
cosine factor by its own derivative in the free value.

The diagonal of every component is **exactly zero**, the diagonal of
\\R\\ being the constant 1. Angle \\\theta\_{ij}\\ belongs to row \\i\\,
so \\\partial_k L\\ is supported on that row alone, though \\\partial_k
R\\ is not.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `s@n_free` symmetric matrices named by `s@free_names`, each
`s@dimension` by `s@dimension` with a zero diagonal.

## See also

[`corr_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_derivative.md),
which assembles it,
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
for the trigonometric derivatives, and
[`param_d2.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.CorrelationParam.md)
for the order above.
