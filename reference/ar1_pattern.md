# The Pattern of an AR(1) Parameter

\\P(\rho)\_{ij} = \rho^{\lvert i-j \rvert}\\ and its four derivatives in
the free value, entry by entry: the composition of a power with the
link.

## Usage

``` r
ar1_pattern(s, sc)
```

## Arguments

- s:

  An
  [`Ar1Param`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- sc:

  The scalars of
  [`econ_scalars`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md).

## Value

A list of five matrices.

## Details

Each distinct lag is composed once and written into every entry that
carries it, the matrix having only \\p\\ distinct values.
