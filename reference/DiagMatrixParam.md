# Diagonal Parameter

The S7 class of diagonal positive matrices, one free value per entry.
Constructed by
[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
or
[`scalar_matrix`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).

## Usage

``` r
DiagMatrixParam(
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

An object of class `DiagMatrixParam`.

## See also

[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scalar_matrix`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)

## Examples

``` r
S7::S7_inherits(diagonal_matrix(3), DiagMatrixParam)
#> [1] TRUE
```
