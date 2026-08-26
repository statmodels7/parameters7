# Log-Determinant Derivatives of an Autoregressive Parameter

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for an
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
parameter, closed form at every order and **separable**: the
log-determinant is a sum with one term per free value, so every mixed
component is exactly zero. At \\q = 2\\ that is 3 of the 6 components at
order 2, 7 of 10 at order 3 and 12 of 15 at order 4, and the zeros are
exact, not merely small.

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector keyed as `param_tuple_names(s, order)` and in that
order: `s@n_free` values at order 1, then `choose(s@n_free + k - 1, k)`
at order \\k\\.

## Details

The four share
[`ar_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_logdet_derivative.md)
and differ only in the order they pass. Compare
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
and
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
whose log-determinants are separable for the same reason, against
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
where a mixed component has content.

## See also

[`ar_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_logdet_derivative.md),
which assembles them, and
[`param_logdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.AutoregressiveParam.md)
for the quantity differentiated.
