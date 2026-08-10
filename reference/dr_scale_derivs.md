# Derivatives of the Inverse Link at Every Scale Coordinate

Returns \\d_j\\ and its first four derivatives in the free value that
carries it, as a matrix with one row per order.

## Usage

``` r
dr_scale_derivs(s, eta)
```

## Arguments

- s:

  A
  [`DrProdParam`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

## Value

A 5 by `p` numeric matrix, rows being orders 0 to 4.
