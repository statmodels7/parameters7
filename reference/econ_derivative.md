# Derivative Components of an Economical Parameter

Assembles one derivative order of \\M = \sigma^2 P(\rho)\\ for a
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
or [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
parameter, from the scale's derivatives and the pattern's. The four
derivative methods of both families are one call each to this function,
differing only in the `pattern` passed.

## Usage

``` r
econ_derivative(s, eta, order, pattern)
```

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  or
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values.

- order:

  The derivative order: 1, 2, 3 or 4.

- pattern:

  A function of the parameter and the scalars of
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md),
  returning a list of five matrices: the pattern \\P(\rho)\\ and its
  four derivatives in the second free value. See
  [`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md)
  and
  [`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md).

## Value

A list of `choose(order + 1, order)`, that is `order + 1`, symmetric
matrices keyed as `param_tuple_names(s, order)` and in that order, each
`s@dimension` by `s@dimension`.

## Details

The value is a **product** of a function of the first free value and a
function of the second, so a component with \\a\\ scale indices and
\\b\\ correlation indices is the \\a\\-th derivative of the scale times
the \\b\\-th derivative of the pattern. Nothing is approximated and no
order is special: only the counts \\a\\ and \\b\\ matter, and a tuple is
fully described by them.

## See also

[`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md)
for the two chains,
[`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md)
for one of the two patterns, and
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
its log-determinant counterpart.
