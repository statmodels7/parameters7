# Log-Determinant Chains of a Correlation Parameter

For each free value, the four derivatives of \\2\log\sin\theta\\ in that
free value: the log-determinant's contribution from one angle.

## Usage

``` r
corr_logdet_chains(s, eta)
```

## Arguments

- s:

  A
  [`CorrelationParam`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A list with one element per free value, each a list of four numbers.

## Details

The factor is triangular, so \\\lvert R \rvert = \prod_i L\_{ii}^2\\ and
\$\$\log\lvert R \rvert = 2 \sum\_{i,k} \log \sin\theta\_{ik},\$\$ a sum
with one term per free value. The log-determinant is therefore
separable, every mixed derivative is exactly zero, and each pure one is
the logarithm composed with the sine table
[`corr_tables`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
already holds.
