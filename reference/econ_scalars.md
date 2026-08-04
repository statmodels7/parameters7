# The Scale and the Correlation of an Economical Parameter

The two scalars a compound-symmetric or AR(1) parameter is built from,
with the first four derivatives of each in its own free value.

## Usage

``` r
econ_scalars(s, eta)
```

## Arguments

- s:

  A
  [`CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  or
  [`Ar1Param`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values.

## Value

A list with `scale` and `rho`, each a list of five numbers: the value
and four derivatives.
