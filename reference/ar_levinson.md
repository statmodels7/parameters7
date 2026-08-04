# The Levinson-Durbin Recursion, Carried in Jets

From the partial autocorrelations as jets, the predictor coefficients of
every order and the autocorrelations out to the dimension, each a jet in
the free values.

## Usage

``` r
ar_levinson(s, r, lay)
```

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- r:

  A list of \\q\\ jets, the partial autocorrelations.

- lay:

  A layout from
  [`jet_layout`](https://statmodels7.github.io/parameters7/reference/jet_layout.md).

## Value

A list with `phi`, the coefficient jets of the last order, and `rho`,
the autocorrelation jets from lag 0 to `dimension - 1`.

## Details

The recursion is sums and products only, so carrying jets through it
gives every derivative to fourth order exactly. The coefficients of the
last order are the autoregression's own, and the autocorrelations beyond
the order come from the Yule-Walker equations.
