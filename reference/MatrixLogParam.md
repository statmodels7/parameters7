# Matrix Logarithm Parameter

The S7 class of unstructured symmetric positive definite matrices in the
matrix logarithm parametrization, \\M = \exp(S)\\ with \\S\\ symmetric
and its lower triangle read straight off the free vector. It is the
other chart onto the same cone
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
parametrizes, with a different set of quantities coming free.

[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
builds one. Every quantity is exact, and every derivative goes through
the eigendecomposition of \\S\\, so the object holds no auxiliary data
at all.

## Usage

``` r
MatrixLogParam(
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

An object of class `MatrixLogParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `positions`, the
row and column of each free value. `n_free` is \\p(p+1)/2\\ and `rank`
is \\p\\.

## See also

[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
the constructor,
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for the other unstructured chart, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# The same cone as log_cholesky(), the same number of free values.
s <- matrix_log(3)
c(matrix_log = s@n_free, log_cholesky = log_cholesky(3)@n_free)
#>   matrix_log log_cholesky 
#>            6            6 

# But different free names: no entry is transformed here.
rbind(matrix_log = s@free_names, log_cholesky = log_cholesky(3)@free_names)
#>              [,1]     [,2]     [,3]     [,4]   [,5]   [,6]  
#> matrix_log   "S1"     "S2"     "S3"     "S2.1" "S3.1" "S3.2"
#> log_cholesky "log_L1" "log_L2" "log_L3" "L2.1" "L3.1" "L3.2"
```
