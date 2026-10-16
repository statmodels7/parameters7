# Log-Determinant of an Autoregressive Parameter

Returns the log-determinant in closed form, from the innovation
variances of the Levinson-Durbin recursion:

\$\$\log\lvert M \rvert = p\log\gamma_0 + \sum\_{k=1}^{q} (p -
k)\log(1 - r_k^{2}).\$\$

The value is a sum of \\q + 1\\ terms at any dimension and is computed
without a factorization or a determinant. Under the rhobit link \\1 -
r_k^2 = \mathrm{sech}^2 z_k\\, and each logarithm is evaluated by
[`log_sech2()`](https://statmodels7.github.io/parameters7/reference/sech2.md)
in the free value \\z_k\\, so the sum keeps its accuracy where \\r_k\\
rounds to \\-1\\ or 1.

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
