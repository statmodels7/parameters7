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

## Value

An object of class `ScaledMatrixParam`.

## See also

[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)

## Examples

``` r
S7::S7_inherits(scaled_matrix(diag(3)), ScaledMatrixParam)
#> [1] TRUE
```
