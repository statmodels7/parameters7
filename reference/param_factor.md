# A Factor of a Parameter's Matrix

Returns the lower triangular \\L\\ with \\M = L L^\top\\, the Cholesky
factor of the matrix the parametrization produces. Simulation is the
usual reason to want it: \\L z\\ with \\z\\ standard normal has
covariance \\M\\, so one factor draws as many vectors as needed. For
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
the factor is what the parametrization holds anyway, so it is returned
without any arithmetic.

## Usage

``` r
param_factor(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  of full rank.

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A `s@dimension` by `s@dimension` lower triangular numeric matrix with a
positive diagonal, satisfying `L %*% t(L) == param_value(s, eta)`.

## Rank deficiency is rejected

A deficient matrix has no Cholesky factor: a triangular \\L\\ with \\L
L^\top = M\\ would need a zero on the diagonal, and the factor is then
not unique. The generic signals an error naming the family and its rank
before dispatching. To simulate from a deficient covariance, take an
eigendecomposition and use the eigenvectors of the non-zero eigenvalues;
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
gives the directions that carry no variance.

## A non-matrix family has no method

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
are checked for and refused by name, their value not being a symmetric
matrix.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns, \\p\\ its side, and \\\eta\\ the free vector.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for a solve through the same factor,
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
for the matrix it factors, and
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
whose free values are the logarithms of this factor's diagonal.

## Examples

``` r
s <- log_cholesky(3)
eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
L <- param_factor(s, eta)
round(L, 4)
#>         [,1]   [,2]   [,3]
#> [1,]  1.1052 0.0000 0.0000
#> [2,]  0.5000 0.8187 0.0000
#> [3,] -0.4000 0.2000 1.3499

# It is a factor, exactly.
max(abs(L %*% t(L) - param_value(s, eta)))
#> [1] 0

# In this parametrization the diagonal is exp() of the first p free values.
all.equal(diag(L), exp(eta[1:3]))
#> [1] TRUE

# What it is for: simulating with the right covariance.
set.seed(1)
y <- t(L %*% matrix(rnorm(3 * 20000), 3, 20000))
round(cov(y) - param_value(s, eta), 2)
#>      v1   v2   v3
#> v1 0.01 0.01 0.00
#> v2 0.01 0.01 0.00
#> v3 0.00 0.00 0.03

# A rank-deficient family has no factor and says so.
try(param_factor(scaled_matrix(crossprod(diff(diag(6), differences = 2))), 0))
#> Error : 'scaled' is rank deficient (4 of 6), so it has no Cholesky factor.
```
