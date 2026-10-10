# Log-Determinant Chains of a Correlation Parameter

Returns, for each free value, the derivative of order `order` of
\\2\log\sin\theta\\ in that free value: the log-determinant's whole
contribution from one angle. The family's four log-determinant
derivative methods read nothing else.

## Usage

``` r
corr_logdet_chains(s, eta, order)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order, an integer from 1 to 4.

## Value

A list of `s@n_free` numbers, the derivative of order `order` of
\\2\log\sin\theta_k\\ in \\\eta_k\\.

## Details

The factor is triangular, so \\\lvert R \rvert = \prod_i L\_{ii}^2\\ and
\$\$\log\lvert R \rvert = 2 \sum\_{i,k} \log \sin\theta\_{ik},\$\$ a sum
with one term per free value. The log-determinant is therefore
separable, every mixed derivative is exactly zero, and each pure one is
the logarithm composed with the sine table that
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
already holds.

The four coefficients \\2(-1)^{j-1}(j-1)!/\sin^j\theta\\ are the
derivatives of \\2\log u\\ at \\u = \sin\theta\\, and
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
chains them onto the sine's own derivatives in the free value, which
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
has already computed.

## See also

[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
for the sine table,
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
for the chain, and
[`corr_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_derivative.md),
the only caller.
