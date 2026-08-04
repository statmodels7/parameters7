# Log-Determinant Components of a Correlation Parameter

Assembles one derivative order of the log-determinant from the per-angle
chains, every mixed component being exactly zero.

## Usage

``` r
corr_logdet_derivative(s, eta, order)
```

## Arguments

- s:

  A
  [`CorrelationParam`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 1 to 4.

## Value

A named numeric vector.
