# Scaled Fixed Matrix Parameter

The S7 class of a fixed symmetric positive semidefinite matrix carried
by a single scale. Constructed by
[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md).

## Usage

``` r
ScaledMatrixParam(
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

  A character vector of length `n_free`, one label per free value. Fixed
  at construction: every consumer builds parameter tables from these.

- param_params:

  A list of whatever the family needs to evaluate itself.

- dimension:

  The side \\p\\ of the matrix.

- rank:

  The rank of the matrix the family produces. Equal to `dimension` for a
  definite family, and strictly less for a rank-deficient precision,
  where it is the dimension of the space the quadratic form penalises.

- null_basis:

  A `dimension` by `dimension - rank` matrix whose columns are an
  orthonormal basis of the null space, or a matrix with no columns when
  the family is of full rank.

- role:

  One of `"covariance"`, `"precision"` or `"either"`. A label: no method
  reads it and no result depends on it. It exists because the name of a
  family does not say which side of a model it parametrises, and the two
  are different models – the inverse of a compound-symmetry matrix is
  compound symmetry, while the inverse of an AR(1) is tridiagonal and
  not AR(1).

## Value

An object of class `ScaledMatrixParam`.

## See also

[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)

## Examples

``` r
S7::S7_inherits(scaled_matrix(diag(3)), ScaledMatrixParam)
#> [1] TRUE
```
