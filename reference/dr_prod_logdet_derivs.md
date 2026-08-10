# Assemble a Scales-Times-Correlation Log-Determinant Derivative

The log-determinant is \\2\sum_j \log d_j + \log\lvert R\rvert\\, so a
component is the scale's own when every index names one scale
coordinate, the correlation's own when every index is a correlation
coordinate, and zero otherwise.

## Usage

``` r
dr_prod_logdet_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`DrProdParam`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order, 1 to 4.

## Value

A named numeric vector keyed as `param_tuple_names(s, order)`.
