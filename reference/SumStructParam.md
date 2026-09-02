# A Non-Negative Combination of Fixed Matrices

The S7 class of a matrix parameter that is a sum of **fixed** symmetric
positive semidefinite matrices, each carried by one positive free value.
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
builds one.

It is the variance-components covariance, and the one family here whose
free values are weights rather than entries: `n_free` is the number of
components, whatever the dimension.

## Usage

``` r
SumStructParam(
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

An object of class `SumStructParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `components`,
`link` and `labels`. `n_free` is the number of components and `rank` is
fixed at construction from their shared null space.

## See also

[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md),
the constructor,
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
for the one-component case with the matrix fixed, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# One free value per component, whatever the side of the matrices.
s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
c(dimension = s@dimension, n_free = s@n_free)
#> dimension    n_free 
#>         3         2 
s@free_names
#> [1] "log_between" "log_within" 
```
