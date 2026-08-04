# Solve of an Autoregressive Parameter

Exact and banded of bandwidth \\q\\. The precision is \\U^\top D^{-1}
U\\ in the prediction form of
[`ar_prediction`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md),
which is what an order-\\q\\ Markov property means: no partial
correlation beyond the lag.

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values.

- b:

  A numeric matrix with `s@dimension` rows.

- ...:

  Unused.

## Value

A numeric matrix.
