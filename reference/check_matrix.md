# Validate a Matrix Handed Back to a Parameter

Checks that `m` is a square symmetric numeric matrix of the parameter's
dimension, and returns it symmetrized. Called by every
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
method on the matrix the caller is asking to invert.

## Usage

``` r
check_matrix(s, m, tol = 1e-08)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, whose `dimension` is read.

- m:

  The matrix supplied by the caller. A non-matrix or non-numeric `m`
  throws, as does a wrong shape (with a message naming the dimension
  required), an `NA` anywhere, and an asymmetry above `tol`.

- tol:

  The relative tolerance for the symmetry check, a single positive
  number. Defaults to `1e-8`, loose enough to accept a matrix assembled
  from a Cholesky factor or a Kronecker product and tight enough to
  reject one that is genuinely not symmetric.

## Value

`m`, symmetrized as `(m + t(m)) / 2`.

## Details

The symmetry check is relative: `m` is rejected when \\\max\_{ij}
\|m\_{ij} - m\_{ji}\|\\ exceeds `tol * max(1, max(abs(m)))`. At the
default a discrepancy of \\10^{-4}\\ throws and one of \\10^{-12}\\
passes. The return value is `(m + t(m)) / 2`, so an asymmetry small
enough to be rounding is averaged away instead of being carried into the
inverse map.
