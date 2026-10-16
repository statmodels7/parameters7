# Log-Determinant of a Correlation Parameter

Computes the log-determinant in closed form. The factor is triangular
with \\L\_{ii} = \prod\_{k\<i} \sin\theta\_{ik}\\ on its diagonal, so
\\\|R\| = \prod_i L\_{ii}^2\\ and

\$\$\log\|R\| = 2 \sum\_{i,k} \log \sin\theta\_{ik},\$\$

one term per free value. It is `2 * sum(log(sines))`, computed without a
factorization, a determinant or an eigendecomposition.

It is always negative or zero, because a correlation matrix has
determinant at most 1, and it equals 0 only at the identity.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number, `0` when `s@n_free` is 0.

## See also

[`param_dlogdet.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CorrelationParam.md)
for its four derivative orders, and
[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md)
for the per-angle chains that they read.
