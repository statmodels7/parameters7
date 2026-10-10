# Log-Determinant Chain of an AR(1) Parameter

Returns a function of the second free value \\z\\ and an order that
gives the order-th derivative of \\q = (p-1)\log(1-\rho^2)\\ in \\z\\,
with \\\rho = \tanh z\\ under the rhobit link.

## Usage

``` r
ar1_logdet_chain(s)
```

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object, whose `dimension` supplies \\p\\.

## Value

A function of `(e, order)`, with `e` a single number and `order` an
integer from 1 to 4, returning a single number.

## Details

\\q = (p-1)\log\mathrm{sech}^2 z\\, so each derivative is \\p - 1\\
times the one that
[`log_sech2_deriv()`](https://statmodels7.github.io/parameters7/reference/sech2.md)
returns, which is written in the free value and keeps its accuracy where
\\\rho\\ rounds to \\-1\\ or 1.

## See also

[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which calls it, and
[`cs_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_chain.md)
for the compound-symmetric counterpart.
