# Unstructured Positive Definite Parameter

The S7 class of unstructured symmetric positive definite matrices in the
log-Cholesky parametrisation. Constructed by
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Usage

``` r
LogCholeskyParam(
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

An object of class `LogCholeskyParam`. Use
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
rather than calling the class directly.

## See also

[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)

## Examples

``` r
S7::S7_inherits(log_cholesky(2), LogCholeskyParam)
#> [1] TRUE
```
