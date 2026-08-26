# The Value a Parameter Produces

Evaluates the map a parametrization stands for: given a free vector
\\\eta\\ on the unconstrained scale, returns the constrained value it
produces. The shape depends on the family. A
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
returns a symmetric \\p \times p\\ matrix,
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
a probability vector, and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
a row-stochastic matrix. Whatever \\\eta\\ is handed in, the value
satisfies the family's constraint, so a caller never has to test the
result.

This is the only method a new family must write. Every derivative order
is then available numerically, and for a matrix family so are the
log-determinant, the solve and the factor.

## Usage

``` r
param_value(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md),
  from any of the constructors:
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
  [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
  and the rest.

- eta:

  A numeric vector of length `s@n_free`, finite in every entry. The
  order of its entries is the order of `s@free_names`. Names are ignored
  and stripped, so a value that arrives labeled by a link does not leak
  that label into the result.

- ...:

  Passed to the method. No method in this package reads it; it is part
  of the signature so that a family written elsewhere can take further
  arguments of its own.

## Value

The constrained value, shaped as the family declares.

- a matrix family:

  a symmetric `s@dimension` by `s@dimension` numeric matrix, with
  `dimnames` `v1`, `v2`, ..., positive definite unless the family
  declares a deficient rank.

- [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md):

  a numeric vector of length `s@n_free + 1`, named `p1`, `p2`, ..., with
  positive entries summing to 1.

- [`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md):

  a square numeric matrix with positive entries and rows summing to 1,
  with `dimnames` `s1`, `s2`, ...

## Validation happens before dispatch

The generic checks `eta` in its own body and passes the checked vector
on, so every method inherits the check, including one written outside
the package. Three things are rejected: a non-numeric `eta`, a length
other than `s@n_free` (the message names both counts and lists the
family's `free_names`), and any `NA`, `NaN` or infinite entry. The last
is a caller error rather than a boundary of the domain, the
unconstrained scale having no edge to reach.

## What is exact and what is not

Every family in this package writes its own `param_value()` out in
closed form. None of them falls back to anything, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
returns `FALSE` for all nine components of all fifteen. The numerical
methods on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
exist for a family written elsewhere.

## Notation

\\\eta\\ is the free vector, the point on the unconstrained scale, of
length \\d = \\ `s@n_free`. \\p\\ is the side of the matrix a matrix
family returns.

## See also

[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
for the inverse map,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
through
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the derivatives of this one,
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
for what a matrix family adds, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
to verify a family end to end.

## Examples

``` r
# An unstructured covariance. Any six numbers give a positive definite
# matrix, which is the property the parametrization exists for.
s <- log_cholesky(3)
M <- param_value(s, c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2))
round(M, 4)
#>         v1      v2      v3
#> v1  1.2214  0.5526 -0.4421
#> v2  0.5526  0.9203 -0.0363
#> v3 -0.4421 -0.0363  2.0221
eigen(M, only.values = TRUE)$values > 0
#> [1] TRUE TRUE TRUE

# Even absurd free values stay inside the cone.
eigen(param_value(s, c(-8, 9, -7, 100, -100, 50)),
      only.values = TRUE)$values > 0
#> [1] TRUE TRUE TRUE

# A probability vector, and a row-stochastic matrix.
sum(param_value(simplex(4), c(1.2, -0.7, 0.3)))
#> [1] 1
rowSums(param_value(transition_matrix(3), rep(0.4, 6)))
#> s1 s2 s3 
#>  1  1  1 

# The three ways a free vector is rejected.
try(param_value(s, c(1, 2)))
#> Error : 'eta' has 2 value(s) but 'log_cholesky' has 6 free value(s) (log_L1, log_L2, log_L3, L2.1, L3.1, L3.2).
try(param_value(s, c(0, 0, 0, 0, 0, Inf)))
#> Error : 'eta' must be finite: the free scale has no boundary to reach.
try(param_value(s, letters[1:6]))
#> Error : 'eta' must be numeric.
```
