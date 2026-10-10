# Cholesky Factorization, With the Rank Decided Before It

Returns the lower triangular Cholesky factor of a symmetric matrix, or
`NULL` when the matrix is not positive definite to the given relative
tolerance. Positive definiteness is decided **before** the factorization
is attempted, from the eigenvalues.

## Usage

``` r
chol_pd(m, tol = 1e-14)
```

## Arguments

- m:

  A symmetric numeric matrix. Not checked for symmetry; the callers pass
  the value of
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  or a matrix already symmetrized by
  [`check_matrix()`](https://statmodels7.github.io/parameters7/reference/check_matrix.md).

- tol:

  The relative tolerance below which the smallest eigenvalue counts as
  zero, so `m` is rejected when `min(ev) <= tol * max(ev)`. Defaults to
  `1e-14`. At that default `diag(c(1, 1e-14))` returns `NULL`, the test
  being an inequality, and `diag(c(1, 1e-13))` factors.

## Value

The lower triangular \\L\\ with \\M = L L^\top\\, or `NULL` when `m` is
empty, has a missing or non-finite entry, has a non-positive largest
eigenvalue, fails the relative test, or makes
[`base::chol()`](https://rdrr.io/r/base/chol.html) signal an error.

## Details

The verdict comes from the spectrum because
[`chol()`](https://rdrr.io/r/base/chol.html) is not a rank test. On a
matrix with an exactly zero eigenvalue the pivot that ought to be zero
comes out positive or negative according to rounding, so
[`chol()`](https://rdrr.io/r/base/chol.html) succeeds on some platforms
and fails on others. The test `min(ev) <= tol * max(ev)` does not depend
on the platform in that way.

The default tolerance, \\10^{-14}\\ (about 45 times the machine
epsilon), is above the rounding level of an eigenvalue that is exactly
zero and below the conditioning that
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
reaches on ordinary free vectors, so a badly conditioned but positive
definite value is accepted.

## See also

[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md),
the generic whose base method calls this, and
[`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)
for the decomposition that the log-determinant uses.
