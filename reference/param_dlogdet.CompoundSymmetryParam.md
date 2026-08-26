# Log-Determinant Derivatives of a Compound Symmetry Parameter

Closed form at all four orders. One page covers
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
because the four differ only in the order taken: \\\log\|M\| =
p\log\sigma^2 + q(\rho)\\ is a **sum** of a function of one free value
and a function of the other, so every mixed component is exactly zero
and the pure ones are the two chains taken separately.

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
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
[`cs_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_terms.md),
chained onto the bounded link.

The four methods return vectors of `order + 1` entries, keyed by the
tuple names of their own order. At \\p = 4\\ and \\\eta = (\log 2,
0.8)\\ the second order is \\(0, -0.856, 0)\\ over
`log_scale:log_scale`, `logit_rho:logit_rho` and `log_scale:logit_rho`:
the first zero from the log link and the last from the separability.

Unlike
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
where the whole log-determinant is linear, the correlation's chain is
non-zero at all four orders, so this family is one of the two where a
check of the higher orders has content.

## See also

[`param_logdet.CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.CompoundSymmetryParam.md)
for the quantity differentiated,
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which assembles all four, and
[`param_dlogdet.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
for the other two-value family.
