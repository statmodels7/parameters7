# Derivatives of a Scales-Times-Correlation Parameter

[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for a
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
parameter. Each component is the scale factor times the correlation's
own component, elementwise, the two groups of free values being
disjoint, so nothing of the correlation family is rederived and nothing
is differenced.

## Arguments

- s:

  A
  [`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

At order 1, a list of `s@n_free` symmetric matrices named by
`s@free_names`; above it, `choose(s@n_free + k - 1, k)` of them keyed as
`param_tuple_names(s, k)` and in that order. Each is `s@dimension` by
`s@dimension` and labeled `v1`, `v2`, ..., `vp` on both margins.

## Details

The four share
[`dr_prod_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_prod_derivs.md)
and differ only in the order they pass. Both factors are sparse, and
where their supports miss each other the component is exactly zero:
measured at \\p = 3\\ that is 1 of the 21 second-order components, 10 of
56 at third order and 37 of 126 at fourth. See
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
for the rule.

## See also

[`dr_prod_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_prod_derivs.md),
which assembles them, and
[`param_dlogdet.DrProdParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DrProdParam.md)
for the log-determinant's own.
