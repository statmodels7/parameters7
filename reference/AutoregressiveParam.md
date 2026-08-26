# Autoregressive Parameter

The S7 class of covariance matrices of \\p\\ consecutive observations of
a stationary autoregression of order \\q\\, parametrized by the marginal
variance and the \\q\\ partial autocorrelations.
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
builds one.

It is the one family here whose free vector does not grow with the
dimension: `n_free` is \\q + 1\\ at \\p = 6\\ and \\q + 1\\ at \\p =
200\\. What grows with \\p\\ is the matrix, not the parametrization.

## Usage

``` r
AutoregressiveParam(
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

- role:

  A single string, one of `"covariance"`, `"precision"` or `"either"`,
  recording which side of a model the matrix parametrizes. **No numeric
  result depends on it.** It is carried because the family name does not
  record it: the same
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  serves either side, and a consumer that prefixes a free name with the
  matrix it describes needs to know which.
  [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
  reads it to give a composite the common role of its blocks, or
  `"either"` when they disagree, and
  [`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
  copies it.

## Value

An object of class `AutoregressiveParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `order`,
`link_scale` and `link_pacf`. `n_free` is \\q + 1\\ and `rank` is
`dimension`.

## See also

[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
the constructor,
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
for the case \\q = 1\\ written out, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# The free vector does not grow with the dimension.
vapply(c(6, 20, 200), function(p) autoregressive(p, order = 2)@n_free,
       numeric(1))
#> [1] 3 3 3

# One scale and q partial autocorrelations, each on its own chart.
autoregressive(20, order = 3)@free_names
#> [1] "log_scale" "z_pacf1"   "z_pacf2"   "z_pacf3"  
```
