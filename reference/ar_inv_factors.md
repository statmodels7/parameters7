# The Prediction Factors of an Inverse Autoregression, With Derivatives

Returns the two factors of \\\Omega = U^\top \mathrm{diag}(\tau) U\\,
each with its derivative arrays to fourth order, stored under the
integer code of the sub-multiset differentiated in.

## Usage

``` r
ar_inv_factors(s, eta, cd)
```

## Arguments

- s:

  An
  [`AutoregressiveInvParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveInvParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`.

- cd:

  The codes of
  [`ar_inv_codes()`](https://statmodels7.github.io/parameters7/reference/ar_inv_codes.md).

## Value

A list with `u` and `tau`, each a list indexed by code: `u` of
`s@dimension` square matrices, `tau` of numeric vectors of that length.

## Details

\\U\\ is unit lower triangular of bandwidth \\q\\, row \\t\\ holding the
coefficients of the best linear predictor of \\y_t\\ from its
predecessors, which for \\t\\ beyond the order are the autoregression's
own. Those lower-order coefficients are the coefficients of the SAME
family at that order, measured exactly:
[`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md)'s
row \\k+1\\ and `autoregressive(p, order = k)`'s `phi` agree to 0. So
the intermediate derivative arrays come from the compiled
Levinson-Durbin recursion run once per order, and nothing is rederived
here. A component differentiating in a partial autocorrelation the order
does not reach is exactly zero.

\\\tau_t = 1/v_t\\ is a PRODUCT of one factor per free value: \\1/v_0\\
from the scale, and \\(1-r_j^2)^{-1}\\ from each partial autocorrelation
the prediction at \\t\\ has reached. A mixed derivative of a product of
univariate factors is the product of their own derivatives, and it is
exactly zero whenever it differentiates in a factor that row does not
carry.

## See also

[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
for the formula and
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
for the recursion the coefficients come from.
