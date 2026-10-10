# Derivative Components of an Autoregressive Parameter

Assembles one derivative order from the kernel of that order: one
Toeplitz matrix per index tuple, filled in compiled code from the column
of autocovariance derivatives that
[`ar_tables()`](https://statmodels7.github.io/parameters7/reference/ar_tables.md)
returns for it. The four methods differ only in the order they pass.

## Usage

``` r
ar_derivative(s, eta, order)
```

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` with dimnames.

## See also

[`ar_tables()`](https://statmodels7.github.io/parameters7/reference/ar_tables.md)
for the recursion, and
[`param_d1.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.AutoregressiveParam.md),
which calls this.
