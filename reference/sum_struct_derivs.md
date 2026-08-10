# Assemble a Sum of Fixed Matrices' Derivatives of a Given Order

The value is linear in the weights, so a component is zero unless every
index of the tuple names the same free value.

## Usage

``` r
sum_struct_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of matrices keyed as `param_tuple_names(s, order)`.
