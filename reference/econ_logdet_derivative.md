# Log-Determinant Components of an Economical Parameter

Assembles one derivative order of \\\log\|M\| = p\log\sigma^2 +
q(\rho)\\ for a
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
or [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
parameter. The four log-determinant derivative methods of both families
are one call each to this function, differing only in the `terms`
passed.

## Usage

``` r
econ_logdet_derivative(s, eta, order, terms)
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

- terms:

  The affine-logarithm terms of \\q\\, as
  [`cs_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_terms.md)
  or
  [`ar1_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_terms.md)
  returns them; see
  [`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md).

## Value

A numeric vector of `order + 1` entries, keyed as
`param_tuple_names(s, order)` and in that order, with every mixed entry
exactly zero.

## Details

The log-determinant is a **sum** of a function of one free value and a
function of the other, where the value is a product, so the structure is
stronger: every mixed component is exactly zero, and the pure ones are
the two chains taken separately. Under the default log link the scale's
own contribution is \\p\eta_1\\, linear, so its second, third and fourth
derivatives vanish too and only the correlation's chain survives above
first order.

## See also

[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
for the correlation's chain,
[`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md)
for the two links' derivatives, and
[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
its counterpart for the matrix itself.
