# Log-Determinant Chains of a Correlation Parameter

Returns, for each free value, the four derivatives of
\\2\log\sin\theta\\ in that free value: the log-determinant's whole
contribution from one angle. The family's four log-determinant
derivative methods read nothing else.

## Usage

``` r
corr_logdet_chains(s, eta)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list of `s@n_free` elements, each a list of four numbers: the first to
fourth derivative of \\2\log\sin\theta_k\\ in \\\eta_k\\.

## Details

The factor is triangular, so \\\lvert R \rvert = \prod_i L\_{ii}^2\\ and
\$\$\log\lvert R \rvert = 2 \sum\_{i,k} \log \sin\theta\_{ik},\$\$ a sum
with one term per free value. The log-determinant is therefore
separable, every mixed derivative is exactly zero, and each pure one is
the logarithm composed with the sine table
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
already holds.

The four coefficients \\2(-1)^{j-1}(j-1)!/\sin^j\theta\\ are the
derivatives of \\2\log u\\ at \\u = \sin\theta\\, and
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
chains them onto the sine's own derivatives in the free value, which
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
has already computed.

## See also

[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
for the sine table,
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
for the chain, and
[`corr_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_derivative.md),
the only caller.
