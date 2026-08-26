# Derivatives of the Inverse Link at Every Scale Coordinate

Returns \\d_j\\ and its first four derivatives in the free value that
carries it, for every \\j\\ at once, as a matrix with one row per order.

## Usage

``` r
dr_scale_derivs(s, eta)
```

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`; only its first \\p\\ entries
  are read.

## Value

A 5 by \\p\\ numeric matrix, row \\k+1\\ holding the \\k\\-th derivative
of the inverse link at each scale coordinate, so row 1 is the scales
themselves.

## Details

Each scale depends on one free value only, so the table is complete:
there are no cross-derivatives between scales to record. It is computed
once per call to
[`dr_prod_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_prod_derivs.md)
and read by
[`dr_scale_factor()`](https://statmodels7.github.io/parameters7/reference/dr_scale_factor.md)
for every component of the order.

The correlation's part of `eta` is ignored, so one table serves every
component whatever the correlation indices are.

## See also

[`dr_scale_factor()`](https://statmodels7.github.io/parameters7/reference/dr_scale_factor.md),
which reads it, and
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md),
the log-determinant's own chain.
