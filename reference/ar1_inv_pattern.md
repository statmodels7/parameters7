# The Pattern of an Inverse AR(1) Parameter

Returns \\G(\rho)\\, the tridiagonal correlation pattern of the
precision of an AR(1), together with its derivatives in the second free
value to the order that the scalars carry.

## Usage

``` r
ar1_inv_pattern(s, z, order)
```

## Arguments

- s:

  An
  [`Ar1InvParam()`](https://statmodels7.github.io/parameters7/reference/Ar1InvParam.md)
  object, whose `dimension` supplies \\p\\.

- z:

  The second free value.

- order:

  The highest derivative order wanted, from 0 to 4.

## Value

A list of `order + 1` square matrices of side `s@dimension`: the pattern
and its derivatives in the second free value.

## Details

The three distinct entries are \\w\\, \\2w-1\\ and \\-\rho w\\ for \\w =
(1-\rho^2)^{-1}\\. Under the rhobit link, with \\\rho = \tanh z\\, they
are \\\cosh^2 z\\, \\\cosh 2z\\ and \\-\tfrac{1}{2}\sinh 2z\\, so one
sequence of derivatives of \\w\\ in \\z\\, from
[`cosh2_derivs()`](https://statmodels7.github.io/parameters7/reference/cosh2_derivs.md),
serves all three: the second is \\2w^{(k)}\\ above order zero and the
third is \\-\tfrac{1}{2}w^{(k+1)}\\. Nothing is chained through the link
and no difference \\1 - \rho^2\\ is formed.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
for the formulas and
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
for the counterpart on the covariance side.
