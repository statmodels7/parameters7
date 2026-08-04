# The Probabilities Behind an Additive Log-Ratio

Declares the probability vector, whose Jacobian is the one of the
softmax map: \\\partial p_i/\partial\eta_j = p_i(\delta\_{ij} - p_j)\\,
the last coordinate contributing \\-p_k p_j\\.

## Arguments

- s:

  A
  [`SimplexParam`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable`](https://statmodels7.github.io/parameters7/reference/param_readable.md).
