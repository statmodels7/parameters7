# The Prediction Form of an Autoregressive Parameter

Returns the unit lower triangular matrix \\U\\ of one-step predictor
coefficients and the vector \\v\\ of innovation variances, which factor
the matrix as \\M = U^{-1} D U^{-\top}\\ with \\D = \mathrm{diag}(v)\\.
It is where
[`param_solve.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_solve.AutoregressiveParam.md)'s
banded precision comes from.

## Usage

``` r
ar_prediction(s, eta)
```

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list with `u`, a unit lower triangular `s@dimension` by `s@dimension`
numeric matrix of bandwidth \\q\\, and `v`, a positive numeric vector of
length `s@dimension`.

## Details

Row \\t\\ holds the coefficients of the best linear predictor of \\y_t\\
from its predecessors, which for \\t\\ beyond the order are the
autoregression's own coefficients, so \\U\\ has bandwidth \\q\\. The
innovation variances fall by a factor \\1 - r_k^2\\ at each of the first
\\q\\ steps and are constant thereafter. At \\p = 6\\, \\q = 2\\,
\\\gamma_0 = 2\\ and partial autocorrelations \\(0.7, -0.3)\\, \\v\\ is
\\(2, 1.02, 0.9282, 0.9282, 0.9282, 0.9282)\\ and every row of \\U\\
from the third on repeats \\(-0.91, 0.30)\\, which are the coefficients
\\\phi = (0.91, -0.30)\\ with a sign.

Both come from the same recursion the derivatives use, run here without
derivatives, since a solve needs no arrays.

## See also

[`param_solve.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_solve.AutoregressiveParam.md),
the only caller, and
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the Levinson-Durbin recursion.
