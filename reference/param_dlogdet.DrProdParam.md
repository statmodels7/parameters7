# Log-Determinant Derivatives of a Scales-Times-Correlation Parameter

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for a
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
parameter. The log-determinant is separable in the scales and separable
from the correlation, so a component mixing two scales, or a scale with
a correlation, is exactly zero: measured at \\p = 3\\, 12 of the 21
second-order components, 43 of 56 at third order and 108 of 126 at
fourth.

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector: at order 1, `s@n_free` values named by `s@free_names`;
above it, `choose(s@n_free + k - 1, k)` values keyed as
`param_tuple_names(s, k)` and in that order.

## Details

The four share
[`dr_prod_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_prod_logdet_derivs.md)
and differ only in the order they pass. At first order the scale entries
are 2 whatever the point, the scales entering \\\log\lvert\Sigma\rvert\\
through \\2\log d_j\\ and the log link canceling its own derivative.

## See also

[`dr_prod_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_prod_logdet_derivs.md),
which assembles them, and
[`param_logdet.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.DrProdParam.md)
for the quantity differentiated.
