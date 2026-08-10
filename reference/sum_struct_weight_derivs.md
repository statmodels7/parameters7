# Derivatives of the Weights of a Sum of Fixed Matrices

Returns each weight and its first four derivatives in the free value
that carries it, as a matrix with one row per order.

## Usage

``` r
sum_struct_weight_derivs(s, eta)
```

## Arguments

- s:

  A
  [`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

## Value

A 5 by `K` numeric matrix, rows being orders 0 to 4.
