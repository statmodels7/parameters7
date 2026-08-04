# Log-Determinant Components of an Autoregressive Parameter

Assembles one derivative order of the log-determinant. It is a sum with
one term per free value, so every mixed component is exactly zero.

## Usage

``` r
ar_logdet_derivative(s, eta, order)
```

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 1 to 4.

## Value

A named numeric vector.
