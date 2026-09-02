# Scaled Fixed Matrix Parameter

The S7 class of a fixed symmetric positive semidefinite matrix \\P\\
carried by a single positive scale, \\M(\eta) = h(\eta) P\\. It is the
only family in the package that is routinely **rank deficient**: \\P\\
may be a difference penalty or a basis Gram matrix with a genuine null
space, and the class records that rank and null basis at construction.

[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
builds one. With `link = NULL` the object holds \\P\\ itself and has no
free value at all, and `param_name` is then `"fixed"` instead of
`"scaled"`.

## Usage

``` r
ScaledMatrixParam(
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

An object of class `ScaledMatrixParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `p`, the fixed
matrix; `link`, the link object or `NULL`; and `logdet_p`, the log
pseudo-determinant of \\P\\ computed once at construction. `n_free` is
1, or 0 when `link` is `NULL`.

## See also

[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
the constructor, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# A ridge is the identity, scaled: full rank, one free value.
r <- scaled_matrix(diag(3))
c(S7::S7_inherits(r, ScaledMatrixParam), rank = r@rank, n_free = r@n_free)
#>          rank n_free 
#>      1      3      1 

# A second-difference penalty is deficient by two, and says so.
q <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
c(dimension = q@dimension, rank = q@rank, null = ncol(q@null_basis))
#> dimension      rank      null 
#>         6         4         2 

# With no link there is nothing to estimate.
f <- scaled_matrix(diag(3), link = NULL)
c(n_free = f@n_free, name = f@param_name)
#>  n_free    name 
#>     "0" "fixed" 
```
