# Log-Determinant Components of an Economical Parameter

Assembles one derivative order of \\\log\lvert M \rvert =
p\log\sigma^2 + q(\rho)\\.

## Usage

``` r
econ_logdet_derivative(s, eta, order, terms)
```

## Arguments

- s:

  The parameter.

- eta:

  A numeric vector of two free values.

- order:

  The derivative order, 1 to 4.

- terms:

  The affine-logarithm terms of \\q\\; see
  [`log_affine_derivs`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md).

## Value

A named numeric vector.

## Details

The log-determinant is a sum of a function of one free value and a
function of the other, so every mixed component is exactly zero and the
pure ones are the two chains taken separately.
