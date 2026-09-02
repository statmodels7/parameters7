# Derivative Arrays of an Inverse AR(1) Parameter

Assembles one order from the product structure \\\Omega = \tau
G(\rho)\\: a component with \\a\\ scale indices and \\b\\ correlation
indices is the \\a\\-th derivative of the reciprocal scale times the
\\b\\-th derivative of the pattern.

## Usage

``` r
ar1_inv_derivative(s, eta, order)
```

## Arguments

- s:

  An
  [`Ar1InvParam()`](https://statmodels7.github.io/parameters7/reference/Ar1InvParam.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A named list of `s@dimension` square matrices, keyed and ordered as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
says.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
for the formulas and
[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
the same assembly on the covariance side.
