# Derivative Components of a Correlation Parameter

Assembles one derivative order by the Leibniz rule on \\R = LL^\top\\.

## Usage

``` r
corr_derivative(s, eta, order)
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

A named list of symmetric matrices.
