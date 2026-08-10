# Assemble a Scales-Times-Correlation Derivative of a Given Order

Multiplies the scale factor by the correlation's own component,
elementwise, for every tuple of the composite's enumeration.

## Usage

``` r
dr_prod_derivs(s, eta, order)
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

A named list of matrices keyed as `param_tuple_names(s, order)`.
