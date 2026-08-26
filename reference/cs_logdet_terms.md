# The Log-Determinant Terms of a Compound-Symmetric Parameter

Returns the affine-logarithm terms of \\q(\rho) = \log\\1 +
(p-1)\rho\\ + (p-1)\log(1-\rho)\\, the correlation's half of the
log-determinant, in the form
[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
consumes. The two terms are the two distinct eigenvalues of the
correlation pattern, the first with multiplicity 1 and the second with
multiplicity \\p-1\\.

## Usage

``` r
cs_logdet_terms(s)
```

## Arguments

- s:

  A
  [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  object, whose `dimension` supplies \\p\\.

## Value

A list of two numeric triples `c(coefficient, a, b)`, standing for
\\c\log(a + b\rho)\\: `c(1, 1, p-1)` and `c(p-1, 1, -1)`.

## See also

[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md),
which differentiates them, and
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which places the result.
