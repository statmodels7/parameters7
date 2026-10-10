# Log-Determinant Chain of a Compound-Symmetric Parameter

Returns a function of the second free value \\\eta_2\\ and an order that
gives the order-th derivative of \\q(\rho) = \log\\1 + (p-1)\rho\\ +
(p-1)\log(1-\rho)\\ in \\\eta_2\\, with the correlation carried by
`bounded_link(-1/(p-1), 1)`.

## Usage

``` r
cs_logdet_chain(s)
```

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object, whose `dimension` supplies \\p\\.

## Value

A function of `(e, order)`, with `e` a single number and `order` an
integer from 1 to 4, returning a single number.

## Details

With \\\omega = (1 + e^{-\eta_2})^{-1}\\ and \\u = 1 - \omega\\, the
correlation is \\\rho = -1/(p-1) + p\omega/(p-1)\\, so \\1 + (p-1)\rho =
p\omega\\ and \\1 - \rho = pu/(p-1)\\ exactly, and

\$\$q = \log p + \log\omega + (p-1)\\\log(p/(p-1)) + \log u\\.\$\$

The two factors are evaluated from \\\eta_2\\ through
`plogis(., log.p = TRUE)`, never from a rounded \\\rho\\, which loses
all accuracy once \\\rho\\ is within rounding of either end of its
interval. The derivatives follow from \\\omega' = \omega u\\ and \\u' =
-\omega u\\: \\u - (p-1)\omega\\, then \\-p\omega u\\, \\-p\omega u(u -
\omega)\\ and \\-p\omega u\\(u - \omega)^2 - 2\omega u\\\\.

## See also

[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which calls it, and
[`ar1_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_chain.md)
for the AR(1) counterpart.
