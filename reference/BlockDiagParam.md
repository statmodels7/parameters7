# Distinct Blocks of a Matrix Parameter

The S7 class of a block-diagonal matrix built from several matrix
parameters, each block a family of its own and each carrying its own
stretch of the free vector.
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
builds one.

It is the composition for independent groups of coefficients whose
structures differ, where
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
is the one for identical blocks sharing a single free vector.

## Usage

``` r
BlockDiagParam(
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

An object of class `BlockDiagParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `blocks`,
`labels`, `rows` and `free` (the ranges each block occupies in the
matrix and in the free vector) and `owner` (the block each free value
belongs to). `dimension`, `n_free` and `rank` are the sums of the
blocks'.

## See also

[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
the constructor,
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
for identical blocks, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# Everything the composite is, is the sum of what the blocks are.
s <- block_diag(subject = log_cholesky(2), time = ar1(3))
c(dimension = s@dimension, n_free = s@n_free, rank = s@rank)
#> dimension    n_free      rank 
#>         5         5         5 

# And the free names carry the block's label, so two blocks of the same
# family stay distinguishable.
block_diag(a = ar1(3), b = ar1(4))@free_names
#> [1] "a_log_scale" "a_z_rho"     "b_log_scale" "b_z_rho"    
```
