# The Levinson-Durbin Recursion With Its Derivatives

Runs the compiled recursion of `ar_taylor_cpp`: the scale and the
partial autocorrelations enter as their link inverses with four
derivatives each, and the autocorrelations, the coefficients and every
partial derivative to fourth order come out as packed arrays.

## Usage

``` r
ar_taylor(s, eta)
```

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A list with `n`, the number of free values; `gamma`, a matrix with one
row per lag; and `phi`, one row per coefficient. Each row packs the
value, then the full derivative tensors of orders one to four, in
row-major order.

## Details

The recursion is sums and products only, so the propagation rules are
the product rule written out per order; every derivative is exact and
nothing is differenced.
