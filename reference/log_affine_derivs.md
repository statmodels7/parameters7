# Derivatives of a Sum of Logarithms of Affine Functions

The first four derivatives at \\r\\ of \\\sum_t c_t \log(a_t + b_t r)\\,
from \\\mathrm{d}^k \log(a + br)/\mathrm{d}r^k = (-1)^{k-1}(k-1)!\\
b^k/(a + br)^k\\.

## Usage

``` r
log_affine_derivs(r, terms)
```

## Arguments

- r:

  The point.

- terms:

  A list of numeric triples `c(coefficient, a, b)`.

## Value

A list of four numbers.

## Details

Both economical families have a log-determinant of this shape – compound
symmetry from its two distinct eigenvalues, AR(1) from \\\lvert R \rvert
= (1-\rho^2)^{p-1}\\ – so the derivatives are written once rather than
transcribed twice.
