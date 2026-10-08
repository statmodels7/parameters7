# The Scale and the Correlation of an Economical Parameter

Returns the two scalars a
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
or [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
parameter is built from, each with its value and its derivatives to
order `order` in its **own** free value. Each scalar depends on one free
value alone, so the two families are separable and their derivative
assembly is a product of two chains. A value alone is `order = 0`, and
no derivative is evaluated then.

## Usage

``` r
econ_scalars(s, eta, order)
```

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  or
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object, whose `param_params$link_scale` and `param_params$link_rho`
  are read.

- eta:

  A numeric vector of two free values.

- order:

  The highest derivative order wanted, an integer from 0 to 4.

## Value

A list with two components, `scale` and `rho`, each a list of
`order + 1` numbers: the value at index 1 and the derivatives in that
free value at the following indices.

## See also

[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md)
and
[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md),
which consume it, and
[`param_value.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_value.CompoundSymmetryParam.md).
