# Identical Blocks of a Matrix Parameter

The S7 class of a block-diagonal matrix built from \\m\\ identical
copies of an inner matrix parameter, \\I_m \otimes S(\eta)\\. The blocks
**share one free vector**, so `n_free` is the inner parameter's however
large \\m\\ is, and the free names are the inner ones unchanged.

It is the first of the four composition wrappers, and the cheapest:
every quantity of the contract is a linear lift of the inner
parameter's, so nothing is rederived and nothing of size \\(md)^2\\ is
decomposed.

## Usage

``` r
KronIdentityParam(
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

An object of class `KronIdentityParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `inner`, the
per-block parameter, and `m`, the number of blocks. `param_name` is
`kron(Im, <inner>)`.

## See also

[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
the constructor,
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
for blocks that are **not** identical and do not share a free vector,
and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# Three copies of a 2 x 2 covariance: a 6 x 6 matrix with three free values.
s <- kron_identity(log_cholesky(2), 3)
S7::S7_inherits(s, KronIdentityParam)
#> [1] TRUE
c(dimension = s@dimension, n_free = s@n_free)
#> dimension    n_free 
#>         6         3 
s@free_names
#> [1] "log_L1" "log_L2" "L2.1"  

# The free vector does not grow with the number of blocks.
vapply(c(1, 5, 100), function(m) kron_identity(log_cholesky(2), m)@n_free,
       integer(1))
#> [1] 3 3 3
```
