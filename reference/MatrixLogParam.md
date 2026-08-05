# Matrix Logarithm Parameter

The S7 class of unstructured symmetric positive definite matrices in the
matrix logarithm parametrization. Constructed by
[`matrix_log`](https://statmodels7.github.io/parameters7/reference/matrix_log.md).

## Usage

``` r
MatrixLogParam(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list(),
  dimension = integer(0),
  rank = integer(0),
  null_basis = integer(0),
  role = character(0)
)
```

## Arguments

- param_name:

  A single character string naming the family.

- n_free:

  The length \\d\\ of the free vector.

- free_names:

  A character vector of length `n_free`.

- param_params:

  A list of whatever the family needs to evaluate itself.

- dimension:

  The side \\p\\ of the matrix.

- rank:

  The rank of the matrix the family produces.

- null_basis:

  A `dimension` by `dimension - rank` matrix whose columns are an
  orthonormal basis of the null space.

- role:

  One of `"covariance"`, `"precision"` or `"either"`. A label: no method
  reads it and no result depends on it.

## Value

An object of class `MatrixLogParam`.

## See also

[`matrix_log`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)

## Examples

``` r
S7::S7_inherits(matrix_log(2), MatrixLogParam)
#> [1] TRUE
```
