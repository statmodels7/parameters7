# Diagonal Parameter

This S7 class represents diagonal matrices with positive entries, and a
linkfunctions7 link carries each entry onto the free scale. Two
constructors return it, and the difference is recorded in
`param_params$shared`:
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
gives one free value per entry, and
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
one free value for the whole diagonal. Every method here reads that
flag, so the two share one class and one set of formulas.

This family combines parameters7 with linkfunctions7: the Jacobian of a
diagonal map is diagonal, which is exactly what a scalar link supplies,
so the link's own exact derivatives to fourth order are used as they
are.

## Usage

``` r
DiagMatrixParam(
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
  the order in which the free vector holds them. Must be unique.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.

- dimension:

  The side \\p\\ of the matrix: a single integer, no smaller than 1. The
  validator rejects a vector and a value below 1.

- rank:

  The rank of the matrix that the family produces, a single integer in
  `0:dimension`. It is a property of the family, so a family whose value
  is positive definite at every \\\eta\\ declares \\p\\ here.

- null_basis:

  A `dimension` by `dimension - rank` numeric matrix whose columns are
  an orthonormal basis of the common null space. Use
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  to obtain one, or `matrix(numeric(0), dimension, 0)` for a full-rank
  family. The validator rejects any other shape.

## Value

An object of class `DiagMatrixParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds two entries:
`link`, the link object, and `shared`, `TRUE` for
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
and `FALSE` for
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md).
`rank` is always \\p\\, because every entry is positive.

## See also

[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
and
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
the two constructors, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties that this class inherits.

## Examples

``` r
# The two constructors return one class, and `shared` distinguishes them.
d <- diagonal_matrix(3)
q <- scalar_matrix(3)
c(S7::S7_inherits(d, DiagMatrixParam), S7::S7_inherits(q, DiagMatrixParam))
#> [1] TRUE TRUE
c(diagonal = d@param_params$shared, scalar = q@param_params$shared)
#> diagonal   scalar 
#>    FALSE     TRUE 
c(diagonal = d@n_free, scalar = q@n_free)
#> diagonal   scalar 
#>        3        1 
```
