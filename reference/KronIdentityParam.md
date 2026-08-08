# Identical Blocks of a Matrix Parameter

The S7 class of a block-diagonal matrix built from \\m\\ identical
copies of an inner matrix parameter. Constructed by
[`kron_identity`](https://statmodels7.github.io/parameters7/reference/kron_identity.md).

## Usage

``` r
KronIdentityParam(
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

  The length \\d\\ of the free vector.

- free_names:

  A character vector of length `n_free`.

- param_params:

  A list of whatever the family needs to evaluate itself.

- dimension:

  The side \\p\\ of the matrix.

- rank:

  The rank of the matrix the family produces.

- null_basis:

  A `dimension` by `dimension - rank` matrix whose columns are an
  orthonormal basis of the null space.

- role:

  One of `"covariance"`, `"precision"` or `"either"`. A label: no method
  reads it and no result depends on it.

## Value

An object of class `KronIdentityParam`.

## See also

[`kron_identity`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)

## Examples

``` r
S7::S7_inherits(kron_identity(log_cholesky(2), 3), KronIdentityParam)
#> [1] TRUE
```
