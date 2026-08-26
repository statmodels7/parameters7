# Solve of an Autoregressive Parameter

Returns \\M^{-1} b\\, exactly and without a factorization: the precision
is \\U^\top D^{-1} U\\ in the prediction form of
[`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md),
assembled and applied.

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- b:

  A numeric matrix with `s@dimension` rows, defaulted to the identity by
  the generic, so `param_solve(s, eta)` is the whole precision.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric matrix with `s@dimension` rows and as many columns as `b`.

## Details

The precision is **banded of bandwidth \\q\\**, which is the content of
an order-\\q\\ Markov property: no partial correlation beyond the lag.
Measured at \\q = 3\\ and \\p = 9\\, all 30 entries outside the band are
exactly 0 and 51 of the 81 entries are non-zero; against
[`solve()`](https://rdrr.io/r/base/solve.html) on the assembled matrix
the agreement is \\7 \times 10^{-16}\\.

The bandedness is not exploited for speed here, the whole \\p\\ by \\p\\
precision being formed and multiplied. What it buys a consumer is the
structure: a model assembling a penalty from this parameter gets a
banded block, and a sparse solver can use it.

## See also

[`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md)
for the factor, and
[`param_logdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.AutoregressiveParam.md),
which reads the same innovation variances.
