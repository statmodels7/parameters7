# The Pattern of an AR(1) Parameter

Returns \\P(\rho)\_{ij} = \rho^{\|i-j\|}\\, the correlation pattern of
an AR(1) matrix, together with its four derivatives in the second free
value. Each entry is a **power** of the correlation, so unlike
[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md)'s
the pattern is not linear and each order needs a genuine chain:
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
composes \\\rho \mapsto \rho^m\\ with the link's own four derivatives.

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
  whose `rho` component supplies the correlation and its four
  derivatives in the free value.

## Value

A list of five `s@dimension` by `s@dimension` matrices: the pattern at
index 1 and its four derivatives at indices 2 to 5, each derivative with
a zero diagonal.

## Details

Only \\p\\ distinct lags occur, so each is composed once and the result
written into every entry that carries it. The lag-0 entries are the
constant 1, so the diagonal of every derivative is exactly zero.

## See also

[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md),
the compound-symmetric counterpart, which is linear in the correlation,
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
for the chain, and
[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
the caller.
