# The Pattern of an AR(1) Parameter

Returns \\P(\rho)\_{ij} = \rho^{\|i-j\|}\\, the correlation pattern of
an AR(1) matrix, together with its derivatives in the second free value
to the order that the scalars carry. Each entry is a **power** of the
correlation, so unlike
[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md)'s
the pattern is not linear and each order needs a genuine chain:
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
composes \\\rho \mapsto \rho^m\\ with the link's own derivatives, one
order at a time.

## Usage

``` r
ar1_pattern(s, sc)
```

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object, whose `dimension` supplies the lags.

- sc:

  The scalars of
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md),
  whose `rho` component supplies the correlation and its derivatives in
  the free value.

## Value

A list of `s@dimension` by `s@dimension` matrices, one more than the
number of derivatives that `sc` carries: the pattern at index 1 and its
derivatives at the following indices, each derivative with a zero
diagonal.

## Details

Only \\p\\ distinct lags occur, so each is composed once and the result
written into every entry that carries it. The lag-0 entries are the
constant 1, so the diagonal of every derivative is exactly zero.

## See also

[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md),
the compound-symmetric counterpart, which is linear in the correlation,
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
for the chain, and
[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
the caller.
