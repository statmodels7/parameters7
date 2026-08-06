# The Scale, the Partial Autocorrelations and the Coefficients

Declares the marginal variance, the partial autocorrelations that
parametrize the family, and the autoregressive coefficients they
produce.

## Arguments

- s:

  An
  [`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable`](https://statmodels7.github.io/parameters7/reference/param_readable.md).

## Details

The coefficients are the quantity an autoregression is usually reported
in and are not obtainable from the covariance the fit prints. They come
out of the Levinson-Durbin recursion the family already propagates
derivatives through, so their derivatives with respect to every free
value are read off the first-order block of those arrays rather than
computed again. The stationary region in the coefficients is not a box,
so their intervals are built on the identity scale: no scalar
transformation expresses the constraint they are under, and an interval
that respected it would not be an interval.
