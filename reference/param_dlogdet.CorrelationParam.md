# Log-Determinant Derivatives of a Correlation Parameter

Closed form at all four orders. One page covers
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for this family because the four differ only in the order of one chain:
\\\log\|R\| = 2\sum\_{i,k}\log\sin\theta\_{ik}\\ is a sum with one term
per free value, so it is **separable**, every mixed component is exactly
zero at every order, and each pure one is the matching derivative of
\\2\log\sin\theta_k\\ in \\\eta_k\\.

[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md)
computes those four derivatives per angle, and
[`corr_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_derivative.md)
places them.

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

A named numeric vector, of the length and keying its own order calls
for; see **Details**.

## Details

The four methods return vectors of different lengths, keyed by the tuple
names of their own order:

- [`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md):
  `s@n_free` entries, keyed by `free_names`.

- [`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md):
  `choose(d + 1, 2)` entries, `param_tuple_names(s, 2)`.

- [`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md):
  `choose(d + 2, 3)` entries, `param_tuple_names(s, 3)`.

- [`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md):
  `choose(d + 3, 4)` entries, `param_tuple_names(s, 4)`.

In each, only the components whose indices are all equal can be
non-zero. At \\p = 3\\ and \\\eta = (0.4, -0.2, 0.6)\\ the second order
is \\(-1.161, -1.215, -1.078)\\ on the three diagonal pairs and exactly
0 on the three off-diagonal ones.

Unlike
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
where the log-determinant is linear and every order above the first
vanishes, this family gives non-zero answers at all four, so it is one
of the families where a check of those orders has content.

## See also

[`param_logdet.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.CorrelationParam.md)
for the quantity differentiated,
[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md)
for the per-angle chains, and
[`param_d3logdet.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.LogCholeskyParam.md)
for a family where these orders vanish.
