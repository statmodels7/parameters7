# Log-Determinant Derivatives of an AR(1) Parameter

Closed form at all four orders. One page covers
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
because the four differ only in the order taken: \\\log\|M\| =
p\log\sigma^2 + (p-1)\log(1-\rho^2)\\ is a **sum** of a function of one
free value and a function of the other, so every mixed component is
exactly zero and the pure ones are the two chains taken separately.

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object.

- eta:

  A numeric vector of two free values, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A named numeric vector of `order + 1` entries, keyed as
`param_tuple_names(s, order)`, with every mixed entry exactly zero.

## Details

The scale's chain is \\p\\ times the derivatives of \\\log h\\, which
under the default log link is \\p\\ at first order and 0 above it. The
correlation's is
[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
on the two terms of
[`ar1_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_terms.md),
chained onto the rhobit link.

The four methods return vectors of `order + 1` entries, keyed by the
tuple names of their own order. At \\p = 4\\ and \\\rho = 0.6\\ the
second order is \\(0, -3.84, 0)\\ over `log_scale:log_scale`,
`z_rho:z_rho` and `log_scale:z_rho`.

Together with
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
this is one of the two families where the higher log-determinant orders
are non-zero, so a check of them has content here where on
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
it would compare two zeros.

## See also

[`param_logdet.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_logdet.Ar1Param.md)
for the quantity differentiated,
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which assembles all four, and
[`param_dlogdet.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
for the sibling family.
