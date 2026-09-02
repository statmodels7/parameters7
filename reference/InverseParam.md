# The Inverse of a Matrix Parameter

The S7 class of the family whose value is the inverse of another
family's.
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
builds one. It carries the inner family's free vector unchanged, so its
coordinates are read as those of the matrix it inverts.

## Usage

``` r
InverseParam(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list(),
  dimension = integer(0),
  rank = integer(0),
  null_basis = integer(0)
)
```

## Arguments

- param_name:

  A single character string naming the family.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer,
  agreeing with `length(free_names)`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order the free vector holds them. Must be unique.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.

- dimension:

  The side \\p\\ of the matrix: a single integer, no smaller than 1. The
  validator rejects a vector and a value below 1.

- rank:

  The rank of the matrix the family produces, a single integer in
  `0:dimension`. It is a property of the family, so a family whose value
  is positive definite at every \\\eta\\ declares \\p\\ here.

- null_basis:

  A `dimension` by `dimension - rank` numeric matrix whose columns are
  an orthonormal basis of the common null space. Use
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  to obtain one, or `matrix(numeric(0), dimension, 0)` for a full-rank
  family. The validator rejects any other shape, and reports both the
  rank and the shape when the two disagree.

## Value

An object of class `InverseParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `inner`, the
family being inverted.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md),
the constructor, and
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
and
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md),
the two families whose inverse has a structure of its own and which are
written out rather than composed.

## Examples

``` r
# The value is the inner family's inverse, and the free vector is the
# inner family's own.
s <- inverse_of(correlation_matrix(3))
s@free_names
#> [1] "z2.1" "z3.1" "z3.2"
eta <- c(0.4, -0.3, 0.8)
max(abs(param_value(s, eta) %*% param_value(correlation_matrix(3), eta) -
        diag(3)))
#> [1] 2.220446e-16
```
