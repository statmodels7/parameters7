# Simplex Parameter

The S7 class of probability vectors on the open simplex, in the additive
log-ratio parametrisation. Constructed by
[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.md).

## Usage

``` r
SimplexParam(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list()
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

An object of class `SimplexParam`.

## See also

[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.md)

## Examples

``` r
S7::S7_inherits(simplex(3), SimplexParam)
#> [1] TRUE
```
