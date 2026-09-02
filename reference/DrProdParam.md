# Scales Times a Correlation Matrix

The S7 class of a covariance written as \\D R D\\, for a diagonal matrix
\\D\\ of positive scales and a correlation matrix \\R\\.
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
builds one.

Its free values are the standard deviations and the correlation's own
coordinates, so the quantities a reader takes off a fitted covariance
are the coordinates themselves.

## Usage

``` r
DrProdParam(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list(),
  dimension = integer(0),
  rank = integer(0),
  null_basis = integer(0)
)
```

## Arguments

- param_name:

  A single character string naming the family.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer,
  agreeing with `length(free_names)`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order the free vector holds them. Must be unique.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.

- dimension:

  The side \\p\\ of the matrix: a single integer, no smaller than 1. The
  validator rejects a vector and a value below 1.

- rank:

  The rank of the matrix the family produces, a single integer in
  `0:dimension`. It is a property of the family, so a family whose value
  is positive definite at every \\\eta\\ declares \\p\\ here.

- null_basis:

  A `dimension` by `dimension - rank` numeric matrix whose columns are
  an orthonormal basis of the common null space. Use
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  to obtain one, or `matrix(numeric(0), dimension, 0)` for a full-rank
  family. The validator rejects any other shape, and reports both the
  rank and the shape when the two disagree.

## Value

An object of class `DrProdParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `cor` (the
correlation parameter), `link` and `p`. `n_free` is \\p\\ plus the
correlation's, and `rank` is \\p\\, this family admitting no deficiency.

## See also

[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md),
the constructor,
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
for the default block, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# The scales come first and the correlation's coordinates follow.
s <- dr_prod(3)
s@free_names
#> [1] "log_sd1" "log_sd2" "log_sd3" "z2.1"    "z3.1"    "z3.2"   

# p standard deviations plus the correlation's p(p-1)/2 angles.
c(p = 4, n_free = dr_prod(4)@n_free)
#>      p n_free 
#>      4     10 
```
