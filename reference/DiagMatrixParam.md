# Diagonal Parameter

The S7 class of diagonal matrices with positive entries, each entry
carried onto the free scale by a linkfunctions7 link. Two constructors
return it, and the difference is recorded in `param_params$shared`:
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
gives one free value per entry, and
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
one free value for the whole diagonal. Every method here reads that
flag, so the two share one class and one set of formulas.

It is the family in which parameters7 and linkfunctions7 meet: the
Jacobian of a diagonal map is diagonal, which is exactly what a scalar
link supplies, so the link's own exact derivatives to fourth order are
used as they are.

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

An object of class `DiagMatrixParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds two entries:
`link`, the link object, and `shared`, `TRUE` for
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
and `FALSE` for
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md).
`rank` is always \\p\\, every entry being positive.

## See also

[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
and
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
the two constructors, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# One class, two constructors, told apart by `shared`.
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
