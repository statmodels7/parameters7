# Transition Matrix Parameter

The S7 class of row-stochastic matrices, each row on the open simplex in
the additive log-ratio parametrisation. Constructed by
[`transition_matrix`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md).

## Usage

``` r
TransitionMatrixParam(
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

  A character vector of length `n_free`.

- param_params:

  A list of whatever the family needs to evaluate itself.

## Value

An object of class `TransitionMatrixParam`.

## See also

[`transition_matrix`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)

## Examples

``` r
S7::S7_inherits(transition_matrix(3), TransitionMatrixParam)
#> [1] TRUE
```
