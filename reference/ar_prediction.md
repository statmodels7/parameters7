# The Prediction Form of an Autoregressive Parameter

The unit lower triangular matrix of one-step predictor coefficients and
the innovation variances, which factor the matrix as \\M = U^{-1} D
U^{-\top}\\.

## Usage

``` r
ar_prediction(s, eta)
```

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A list with the matrix `u` and the vector `v`.

## Details

Row \\t\\ holds the coefficients of the best linear predictor of \\y_t\\
from its predecessors, which for \\t\\ beyond the order are the
autoregression's own coefficients, so \\U\\ has bandwidth \\q\\. The
innovation variances fall by a factor \\1 - r_k^2\\ at each of the first
\\q\\ steps and are constant thereafter.
