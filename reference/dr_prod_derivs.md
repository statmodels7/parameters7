# Assemble a Scales-Times-Correlation Derivative of a Given Order

Splits each tuple of the composite's enumeration into its scale indices
and its correlation indices, and multiplies
[`dr_scale_factor()`](https://statmodels7.github.io/parameters7/reference/dr_scale_factor.md)
by the correlation's own component of the matching order, elementwise.

## Usage

``` r
dr_prod_derivs(s, eta, order)
```

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.

## Details

The correlation's arrays are fetched at most once per order and re-keyed
by the sorted local index tuple, the same device
[`block_derivs_by_tuple()`](https://statmodels7.github.io/parameters7/reference/block_derivs_by_tuple.md)
uses and for the same reason: the composite's names are not the block's,
but the sorted index tuple is a key both sides can compute.

A component whose correlation part is empty reads the correlation's
**value**, which is why
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
appears in a derivative routine.

## See also

[`dr_scale_factor()`](https://statmodels7.github.io/parameters7/reference/dr_scale_factor.md)
for one factor, and
[`param_d1.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.DrProdParam.md),
which calls this.
