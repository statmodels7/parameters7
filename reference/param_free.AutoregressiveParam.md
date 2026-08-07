# Free Vector of an Autoregressive Parameter

The marginal variance is the common diagonal entry and the partial
autocorrelations come from the Levinson-Durbin recursion run forwards on
the autocorrelations. The matrix is then checked against the pattern
those values imply, so a Toeplitz matrix that is not the covariance of
an autoregression of this order is rejected rather than fitted.

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- m:

  A covariance matrix of a stationary autoregression.

- ...:

  Unused.

## Value

A named numeric vector of free values.
