# The Precision of an Autoregression of Order q

The S7 class of the family whose value is the **inverse** of an
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
matrix: banded of bandwidth \\q\\, the precision of an autoregression of
that order.
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
builds one.

## Usage

``` r
AutoregressiveInvParam(
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

An object of class `AutoregressiveInvParam`, a subclass of
[`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
adding no properties of its own. `param_params` holds `inner`, the
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
family being inverted.

## See also

[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md),
the constructor, and
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
for the order-one case.

## Examples

``` r
s <- autoregressive_inv(6, order = 2)
s@free_names
#> [1] "log_scale" "z_pacf1"   "z_pacf2"  
```
