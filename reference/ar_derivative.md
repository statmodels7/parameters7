# Derivative Components of an Autoregressive Parameter

Assembles one derivative order by reading the matching component out of
every jet.

## Usage

``` r
ar_derivative(s, eta, order)
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

A named list of symmetric matrices.
