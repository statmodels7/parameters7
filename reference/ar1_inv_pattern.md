# The Pattern of an Inverse AR(1) Parameter

Returns \\G(\rho)\\, the tridiagonal correlation pattern of the
precision of an AR(1), together with its four derivatives in the second
free value.

## Usage

``` r
ar1_inv_pattern(s, sc)
```

## Arguments

- s:

  An
  [`Ar1InvParam()`](https://statmodels7.github.io/parameters7/reference/Ar1InvParam.md)
  object, whose `dimension` supplies \\p\\.

- sc:

  The scalars of
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md)
  read on the inner
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
  whose `rho` entry supplies the correlation and its link's derivatives.

## Value

A list of five `s@dimension` square matrices: the pattern and its four
derivatives in the second free value.

## Details

The three distinct entries are \\w\\, \\2w-1\\ and \\-\rho w\\ for \\w =
(1-\rho^2)^{-1}\\, so one sequence of derivatives of \\w\\ in \\\rho\\
serves all three. Writing \\w\\ by partial fractions as
\\\tfrac{1}{2}\\(1-\rho)^{-1} + (1+\rho)^{-1}\\\\ makes every order an
exact expression rather than a repeated quotient rule, and it is finite
throughout \\\|\rho\| \< 1\\, which the rhobit link guarantees.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
for the formulas and
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
for the counterpart on the covariance side.
