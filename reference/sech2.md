# The Squared Hyperbolic Secant and Its Logarithm

`sech2()` returns \\\mathrm{sech}^2 z = 1 - \tanh^2 z\\, `log_sech2()`
its logarithm, and `log_sech2_deriv()` the derivative of the logarithm
of order 1 to 4 in \\z\\. Under the rhobit link these are the factor
\\1 - \rho^2\\ of
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and the factors \\1 - r_k^2\\ of
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
written in the free value.

## Usage

``` r
sech2(z)

log_sech2(z)

log_sech2_deriv(z, order)
```

## Arguments

- z:

  The free value, a numeric vector.

- order:

  The derivative order, an integer from 1 to 4.

## Value

A numeric vector of the length of `z`.

## Details

\\\mathrm{sech}^2 z\\ is evaluated as \\4a/(1 + a)^2\\ with \\a =
e^{-2\|z\|}\\, and its logarithm as \\2\log 2 - 2\|z\| - 2\log(1 + a)\\.
Neither is computed from \\1 - \tanh^2 z\\, which loses all accuracy
once \\\tanh z\\ rounds to \\-1\\ or 1. With \\t = \tanh z\\ and \\u =
\mathrm{sech}^2 z\\, from \\t' = u\\ and \\u' = -2tu\\, the four
derivatives of the logarithm are \\-2t\\, \\-2u\\, \\4tu\\ and \\4u(u -
2t^2)\\.

## See also

[`ar1_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_chain.md)
and
[`ar_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_logdet_derivative.md),
which call `log_sech2_deriv()`.
