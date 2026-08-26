# Solve Through a Parameter's Matrix

Returns \\M^{-1} B\\ for a right-hand side \\B\\, computed through a
factorization instead of by forming the inverse. This is what a Gaussian
quadratic form needs: `crossprod(r, param_solve(s, eta, r))` is \\r^\top
M^{-1} r\\ at the cost of one triangular solve, where inverting and
multiplying would cost more and be less accurate. Called with no `B` it
does return the inverse, which is convenient for reading a covariance
off a precision.

## Usage

``` r
param_solve(s, eta, b = NULL, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  of full rank.

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- b:

  A numeric matrix or vector with `s@dimension` rows; a vector is
  treated as a one-column matrix. Defaults to `NULL`, which stands for
  the identity and returns the inverse. A wrong number of rows throws a
  message naming the number required.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A numeric matrix with `s@dimension` rows and as many columns as `b`, so
`s@dimension` by `s@dimension` when `b` is left out. A vector `b`
returns a one-column matrix, not a vector.

## Rank deficiency is rejected

A deficient family signals an error and names its rank, instead of
returning a pseudo-inverse. What a consumer of an improper prior needs
is the quadratic form and the log pseudo-determinant: penalized normal
equations invert \\X^\top X + \lambda P\\, which is non-singular even
where \\P\\ is not, and the consumer assembles that matrix itself. A
pseudo-inverse returned here would be a plausible matrix answering a
question nobody asked, and the caller would have no way to tell.

## A non-matrix family has no method

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
inherit
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
alone, not
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
and the generic checks for that before dispatching, so the message names
the family and says its value is not a symmetric matrix.

## Cost

The base method on
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
which most families take, works through a Cholesky factor: \\O(p^3)\\
once plus \\O(p^2)\\ per column of \\B\\. Seven families override it
with a closed-form inverse.
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
writes its tridiagonal precision out entry by entry,
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
uses Sherman-Morrison, and
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
and
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
solve blockwise, so none of them decomposes anything of side \\p\\.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns, \\p\\ its side, and \\\eta\\ the free vector.

## See also

[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for the factor itself,
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
for the other half of a Gaussian log-density, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for what a deficient family offers in place of an inverse.

## Examples

``` r
s <- log_cholesky(3)
eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)

# With no b it is the inverse, and agrees with solve().
max(abs(param_solve(s, eta) - solve(param_value(s, eta))))
#> [1] 4.440892e-16

# The quadratic form of a Gaussian log-density, in one solve.
r <- c(0.4, -1.1, 0.7)
q <- drop(crossprod(r, param_solve(s, eta, r)))
c(q, drop(r %*% solve(param_value(s, eta)) %*% r))
#> [1] 3.314441 3.314441

# A vector b comes back as a one-column matrix.
dim(param_solve(s, eta, r))
#> [1] 3 1

# A rank-deficient family refuses, and says what to do instead.
r_def <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
try(param_solve(r_def, 0))
#> Error : 'scaled' is rank deficient (4 of 6), so it has no inverse. A consumer of
#>   an improper prior needs the quadratic form and the log
#>   pseudo-determinant, not a pseudo-inverse: assemble the matrix the
#>   model actually inverts and solve that.

# So does a family whose value is not a symmetric matrix.
try(param_solve(simplex(3), c(0, 0)))
#> Error : param_solve() is not defined for 'simplex': the value is not a symmetric matrix.
```
