# The Log-Determinant Terms of an AR(1) Parameter

Returns the affine-logarithm terms of \\q(\rho) = (p-1)\log(1-\rho^2)\\,
the correlation's half of the log-determinant, in the form
[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
consumes. The quadratic is split as \\(p-1)\\\log(1-\rho) +
\log(1+\rho)\\\\ so that both pieces are logarithms of functions
**affine** in the correlation, which is the only shape
[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
differentiates.

## Usage

``` r
ar1_logdet_terms(s)
```

## Arguments

- s:

  An
  [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  object, whose `dimension` supplies \\p\\.

## Value

A list of two numeric triples `c(coefficient, a, b)`, standing for
\\c\log(a + b\rho)\\: `c(p-1, 1, -1)` and `c(p-1, 1, 1)`.

## See also

[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md),
which differentiates them,
[`cs_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_terms.md)
for the compound-symmetric pair, and
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which places the result.
