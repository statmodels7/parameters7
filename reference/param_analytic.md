# Which Derivatives of a Parameter Are Analytic

Returns, for each derivative generic that
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports, `TRUE` where the quantity is computed in closed form. It
differs from `!param_is_numerical(s)` in two ways. The base-class
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
and
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
are identities in the derivative arrays, so they count as analytic where
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
for the second, are. And a composite family, which registers methods of
its own that are built from the parameters held in its `param_params`,
counts a quantity as analytic only where every such parameter's quantity
is analytic too.

## Usage

``` r
param_analytic(s)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

## Value

A named logical vector with the names of
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md).

## See also

[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md),
which declares the analytic quantities to the stencil, and
[`param_children()`](https://statmodels7.github.io/parameters7/reference/param_children.md).
