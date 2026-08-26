# Cholesky Factorization, With the Rank Decided Before It

Returns the lower triangular Cholesky factor of a symmetric matrix, or
`NULL` when the matrix is not positive definite to the given relative
tolerance. Positive definiteness is decided **before** the factorization
is attempted, from the eigenvalues.

## Usage

``` r
chol_pd(m, tol = 1e-12)
```

## Arguments

- m:

  A symmetric numeric matrix. Not checked for symmetry; the callers have
  already symmetrized.

- tol:

  The relative tolerance below which the smallest eigenvalue counts as
  zero, so `m` is rejected when `min(ev) <= tol * max(ev)`. Defaults to
  `1e-12`. At that default `diag(c(1, 1e-12))` returns `NULL`, the test
  being an inequality, and `diag(c(1, 1e-11))` factors.

## Value

The lower triangular \\L\\ with \\M = L L^\top\\, or `NULL` when `m` has
no eigenvalues, has an `NA` among them, has a non-positive largest
eigenvalue, fails the relative test, or makes
[`base::chol()`](https://rdrr.io/r/base/chol.html) raise after all.

## Details

The verdict comes from the spectrum, because
[`chol()`](https://rdrr.io/r/base/chol.html) is not a rank test. On a
matrix with an exactly zero eigenvalue the pivot that ought to be zero
comes out positive or negative according to rounding, so
[`chol()`](https://rdrr.io/r/base/chol.html) succeeds on some platforms
and fails on others, and a branch that asks it whether a matrix is
usable gets a different answer on different machines.
`min(ev) <= tol * max(ev)` is a statement about the matrix; a caught
error is a statement about the arithmetic.

## See also

[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md),
the generic whose base method calls this, and
[`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)
for the decomposition the deficient branches use instead.
