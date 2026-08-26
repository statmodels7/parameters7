# Log-Determinant of an Autoregressive Parameter

Closed form, from the innovation variances of the Levinson-Durbin
recursion:

\$\$\log\lvert M \rvert = p\log\gamma_0 + \sum\_{k=1}^{q} (p -
k)\log(1 - r_k^{2}).\$\$

A sum of \\q + 1\\ terms at any dimension, with no factorization and no
determinant taken. Measured against
[`determinant()`](https://rdrr.io/r/base/det.html), the gap is \\2
\times 10^{-15}\\ at \\p = 5\\ and \\4 \times 10^{-14}\\ at \\p = 100\\.

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

A single number.

## See also

[`param_dlogdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
for its derivatives, and
[`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md),
which returns the same innovation variances.
