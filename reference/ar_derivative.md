# Derivative Components of an Autoregressive Parameter

Assembles one derivative order by running
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
once and reading the matching column out of the packed arrays for each
tuple of the order. The four methods differ only in the order they pass.

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

[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
for the recursion,
[`ar_assemble()`](https://statmodels7.github.io/parameters7/reference/ar_assemble.md)
for one component, and
[`param_d1.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.AutoregressiveParam.md),
which calls this.
