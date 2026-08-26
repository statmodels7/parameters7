# The Pattern of a Compound-Symmetric Parameter

Returns \\P(\rho) = I + \rho(J - I)\\, the correlation pattern of a
compound-symmetric matrix, together with its four derivatives in the
second free value. The pattern is **linear** in the correlation, so
every derivative is the matching derivative of \\\rho\\ times the
constant matrix \\J - I\\, and no order needs its own algebra.

## Usage

``` r
cs_pattern(s, sc)
```

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object, whose `dimension` supplies \\I\\ and \\J\\.

- sc:

  The scalars of
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md),
  whose `rho` component supplies the correlation and its four
  derivatives.

## Value

A list of five `s@dimension` by `s@dimension` matrices: the pattern at
index 1 and its four derivatives at indices 2 to 5. Each derivative has
a zero diagonal, the diagonal of the pattern being the constant 1.

## See also

[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md),
the caller, and
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md),
the AR(1) counterpart, which is not linear in the correlation.
