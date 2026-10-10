# The Prediction Factors of an Inverse Autoregression, With Derivatives

Returns the two factors of \\\Omega = U^\top \mathrm{diag}(\tau) U\\,
each with its derivative arrays to fourth order, stored under the
integer code of the sub-multiset differentiated in.

## Usage

``` r
ar_inv_factors(s, eta, cd, order)
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

- order:

  The derivative order the codes belong to, which bounds every factor's
  order.

## Value

A list with `u` and `tau`, each a list indexed by code: `u` of
`s@dimension` square matrices, `tau` of numeric vectors of that length.

## Details

\\U\\ is unit lower triangular of bandwidth \\q\\, row \\t\\ holding the
coefficients of the best linear predictor of \\y_t\\ from its
predecessors, which for \\t\\ beyond the order are the autoregression's
own. Those lower-order coefficients are the coefficients of the same
family at that order (row \\k+1\\ of
[`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md),
negated and read backwards, is the `phi` of
`autoregressive(p, order = k)` up to rounding), so the intermediate
derivative arrays come from the compiled Levinson-Durbin recursion run
once per order, and nothing is rederived here. A component
differentiating in a partial autocorrelation that the order does not
reach is exactly zero.

\\\tau_t = 1/v_t\\ is a product of one factor per free value: \\1/v_0\\
from the scale, and \\(1-r_j^2)^{-1} = \cosh^2 z_j\\ from each partial
autocorrelation that the prediction at \\t\\ has reached, differentiated
in \\z_j\\ by
[`cosh2_derivs()`](https://statmodels7.github.io/parameters7/reference/cosh2_derivs.md).
A mixed derivative of a product of univariate factors is the product of
their own derivatives, and it is exactly zero whenever it differentiates
in a factor that row does not carry.

## See also

[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
for the formula and
[`ar_tables()`](https://statmodels7.github.io/parameters7/reference/ar_tables.md)
for the recursion the coefficients come from.
