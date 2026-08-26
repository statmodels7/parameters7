# The Free Vector Behind a Value

Inverts the map: given a value in the family's set, returns the free
vector \\\eta\\ that
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
would send there. The two are a bijection onto the set, so the round
trip closes to machine precision, and every one of the fifteen families
in this package inverts exactly, with a worst measured error of \\5
\times 10^{-16}\\.

Use it to start an optimizer from a matrix rather than from a free
vector: fit an unstructured covariance by moments, invert it, and hand
the result in as a starting point.

## Usage

``` r
param_free(s, m, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- m:

  A value of the family's shape: a symmetric `s@dimension` by
  `s@dimension` numeric matrix on the matrix branch, a probability
  vector for
  [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md),
  a row-stochastic matrix for
  [`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md).
  It must satisfy the family's constraint; see **Details** for what each
  family rejects.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`, such
that `param_value(s, param_free(s, m))` recovers `m`.

## Exact or refused, never approximated

A family either writes its inverse out or signals an error. The base
method on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
does the second on behalf of a family that has not written the first,
naming the family in the message. An inverse found by minimizing
\\\lVert V(\eta) - m \rVert\\ would return a plausible \\\eta\\ for a
matrix that is not in the set at all, and the caller could not tell that
answer from a correct one.

## A value outside the set is rejected

Every method checks `m` against the family's own constraint before
inverting, and the message says what failed.
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
rejects a matrix that is not positive definite, with the verdict taken
from the eigenvalues;
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
rejects a vector that does not sum to one, and does not renormalize it,
because a silent repair would hide the caller's mistake. The shared
checks are the shape and the symmetry: a matrix of the wrong side is
rejected with a message naming the side required, and an asymmetry above
\\10^{-8}\\ relative is rejected while one below it is averaged away.

## Notation

\\\eta\\ is the free vector, of length \\d = \\ `s@n_free`, and \\m\\
the value on the constrained scale.

## See also

[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
the map this inverts, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which closes the round trip as one of its checks.

## Examples

``` r
# The round trip closes exactly, in both directions.
s <- log_cholesky(3)
eta <- c(0.3, -0.2, 0.5, 0.1, -0.4, 0.2)
param_free(s, param_value(s, eta))
#> log_L1 log_L2 log_L3   L2.1   L3.1   L3.2 
#>    0.3   -0.2    0.5    0.1   -0.4    0.2 
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 5.551115e-17

# Starting an optimizer from a matrix: invert the sample covariance.
set.seed(1)
y <- matrix(rnorm(300), 100, 3) %*% chol(matrix(c(1, .5, .2, .5, 1, .3,
                                                  .2, .3, 1), 3, 3))
start <- param_free(s, cov(y))
round(start, 3)
#> log_L1 log_L2 log_L3   L2.1   L3.1   L3.2 
#> -0.107 -0.187 -0.017  0.448  0.198  0.172 
max(abs(param_value(s, start) - cov(y)))
#> [1] 1.110223e-16

# A matrix outside the set is rejected, and the message says why.
try(param_free(s, diag(c(1, 1, -1))))
#> Error : 'm' is not positive definite, so it is not in the set log_cholesky()
#>   parametrizes. The verdict is spectral, not a failed factorization.

# So is a vector off the simplex; it is not renormalized.
try(param_free(simplex(3), c(0.5, 0.6, 0.2)))
#> Error : 'm' does not sum to one, so it is not on the simplex. It is rejected
#>   rather than renormalized, because a silent repair would mask the
#>   caller's defect.
```
