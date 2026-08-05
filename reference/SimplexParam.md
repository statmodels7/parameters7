# Simplex Parameter

The S7 class of probability vectors on the open simplex, in the additive
log-ratio parametrization. Constructed by
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
  at construction: every consumer builds parameter tables from these. A
  label names the coordinate rather than the quantity the coordinate
  produces, and the families here follow one convention for it. Where a
  link carries a constrained quantity onto the free scale, the label
  records that link, so a variance appears as `"log_scale"` and a
  correlation as `"z_rho"`; where the coordinate is already unrestricted
  the label is the plain name of the quantity, as the below-diagonal
  entries `"L2.1"` of a Cholesky factor are. The distinction matters
  outside the family: a consumer flattens the free vector into scalar
  parameters carrying identity links, so a label promising a bounded
  quantity reports a number on a scale that number is not on.

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
