# The Column of a Packed Derivative Record

Where a derivative component sits in a row of
[`ar_taylor`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)'s
output: the value first, then the tensors of orders one to four in
row-major order.

## Usage

``` r
ar_pack_col(n, order = 0L, tuple = NULL)
```

## Arguments

- n:

  The number of free values.

- order:

  The derivative order, or 0 for the value.

- tuple:

  The index tuple, 1-based.

## Value

A single column index.
