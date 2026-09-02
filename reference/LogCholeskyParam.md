# Unstructured Positive Definite Parameter

The S7 class of unstructured symmetric positive definite matrices in the
log-Cholesky parametrization, \\M = L L^\top\\ with \\L\\ lower
triangular and positive on the diagonal. It is the class the methods of
that family dispatch on: given only
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
everything else here is closed form, including all four derivative
orders and all four orders of the log-determinant.

Call
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
to build one. The class constructor takes the eight properties directly
and does no work; using it means computing `free_names`, `rank` and
`param_params` by hand.

## Usage

``` r
LogCholeskyParam(
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

An object of class `LogCholeskyParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. It carries `dimension`, `rank` (always
\\p\\), `null_basis` (\\p\\ by 0), `param_name` (`"log_cholesky"`),
`n_free` (\\p(p+1)/2\\), `free_names` and `param_params`, whose only
entry is `positions`.

## See also

[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
the constructor, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# The class is what dispatch keys on.
s <- log_cholesky(2)
c(S7::S7_inherits(s, LogCholeskyParam),
  S7::S7_inherits(s, matrix_parameter))
#> [1] TRUE TRUE

# Belonging to it is what gives the family its closed forms: nothing here
# is obtained by finite differences.
any(param_is_numerical(s))
#> [1] FALSE
```
