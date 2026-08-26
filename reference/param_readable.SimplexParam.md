# The Probabilities Behind an Additive Log-Ratio

Declares the probability vector itself, `p1` ... `pK`, which a reader of
a mixture weight or a state distribution needs: the free values are log
ratios against the last category and are not interpretable on their own.

The Jacobian is the softmax map's, \\\partial p_i/\partial\eta_j =
p_i(\delta\_{ij} - p_j)\\, with the reference category's row
contributing \\-p_K p_j\\. It is \\K\\ by \\K-1\\: one more quantity
than there are free values, since the reference probability is
determined by the others. Its **columns** sum to zero, \\\sum_i p_i\\
being the constant 1, measured at \\7 \times 10^{-18}\\.

Every interval is built on the logit scale, so it stays inside \\(0,
1)\\.

## Arguments

- s:

  A
  [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  object.

- eta:

  A numeric vector of free values, of length \\K-1\\.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md).
