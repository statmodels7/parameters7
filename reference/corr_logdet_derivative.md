# Log-Determinant Components of a Correlation Parameter

Assembles one derivative order of the log-determinant from the per-angle
chains of
[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md).
Every mixed component is exactly zero, the log-determinant being a sum
with one term per free value, and each pure one is that angle's chain at
the matching order. The four log-determinant derivative methods of the
family are one call each to this function.

## Usage

``` r
corr_logdet_derivative(s, eta, order)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A numeric vector of `choose(s@n_free + order - 1, order)` entries, keyed
as `param_tuple_names(s, order)` and in that order, or an empty named
vector when there are no free values.

## See also

[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md)
for the chains, and
[`param_dlogdet.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CorrelationParam.md),
which documents all four orders.
