# Compound Symmetry Parameter

The S7 class of compound-symmetric covariance matrices: equal variances
and one common correlation. Constructed by
[`compound_symmetry`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md).

## Usage

``` r
CompoundSymmetryParam(
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

An object of class `CompoundSymmetryParam`.

## See also

[`compound_symmetry`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)

## Examples

``` r
S7::S7_inherits(compound_symmetry(3), CompoundSymmetryParam)
#> [1] TRUE
```
