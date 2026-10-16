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

The precision is banded of bandwidth \\q\\, which expresses an
order-\\q\\ Markov property: partial correlations vanish beyond lag
\\q\\. Every entry outside the band is exactly 0.

The bandedness is not exploited for speed here, because the whole \\p\\
by \\p\\ precision is formed and multiplied. The bandedness is still
useful to a consumer as structure: a model that assembles a penalty from
this parameter receives a banded block, which a sparse solver can use.

## See also

[`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md)
for the factor, and
[`param_logdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_logdet.AutoregressiveParam.md),
which reads the same innovation variances.
