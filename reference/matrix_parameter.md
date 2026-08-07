# Constrained Symmetric Matrix Parameter

The abstract S7 class of the symmetric positive semidefinite branch: a
[`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
whose value is a symmetric matrix, together with the quantities only a
matrix can answer – the rank, the null space, the
log-(pseudo-)determinant, the solve and the factor.

## Usage

``` r
matrix_parameter(
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

An object of class `matrix_parameter`. The class is abstract; use one of
the constructors, such as
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Details

The rank and the null space are properties of the family rather than of
a point: scaling a matrix by a positive number and summing positive
semidefinite matrices both leave the null space alone. They are computed
once, at construction, from the components rather than from an assembled
matrix, because the numerical determination of a rank from an assembled
matrix is not scale invariant while the null space is. See
[`param_null_basis`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md).

A parameter that is not a matrix –
[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.md),
a
[`transition_matrix`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
– inherits from
[`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
directly, so
[`param_logdet`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
and
[`param_solve`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
do not exist for it by construction rather than by a run-time rejection.

## See also

[`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md),
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)

## Examples

``` r
S7::S7_inherits(log_cholesky(3), matrix_parameter)
#> [1] TRUE
```
