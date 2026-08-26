# The Levinson-Durbin Recursion With Its Derivatives

Runs the compiled recursion `ar_taylor_cpp`. The scale and the partial
autocorrelations enter as their link inverses carrying four derivatives
each, and the autocorrelations, the coefficients and every partial
derivative to fourth order come back as packed arrays.

## Usage

``` r
ar_taylor(s, eta)
```

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list with `n`, the number of free values \\q + 1\\; `gamma`, a matrix
with one row per lag \\0, \dots, p-1\\; and `phi`, one row per
autoregressive coefficient. Each row packs the value first and then the
full derivative tensors of orders one to four in row-major order, so it
has \\1 + n + n^2 + n^3 + n^4\\ columns.

## Details

The recursion is sums and products only, so the propagation rules are
the product rule written out per order: every derivative is exact and
nothing is differenced.

It always fills **all four orders**, whatever order the caller went on
to want. Measured at \\q = 2\\, `gamma` comes back 6 by 121 for \\p =
6\\, and \\121 = 1 + 3 + 3^2 + 3^3 + 3^4\\. That is why
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
costs nearly what
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
costs (0.00516 s against 0.01203 s at \\q = 2\\, \\p = 200\\): the
difference between them is the R-level assembly, not the recursion.

## See also

[`ar_pack_col()`](https://statmodels7.github.io/parameters7/reference/ar_pack_col.md)
for the column a component sits in,
[`ar_assemble()`](https://statmodels7.github.io/parameters7/reference/ar_assemble.md)
for the matrix built from one column, and
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the recursion itself.
