# Constrained Matrix Parameter

The abstract S7 class of covariance structures: a map from an
unconstrained vector \\\eta \in \mathbb{R}^{d}\\ to a symmetric matrix
in some constrained set, together with its derivatives and the
quantities a likelihood asks of the matrix.

## Usage

``` r
covstruct(
  struct_name = character(0),
  dimension = integer(0),
  n_free = integer(0),
  free_names = character(0),
  rank = integer(0),
  null_basis = integer(0),
  role = character(0),
  struct_params = list()
)
```

## Arguments

- struct_name:

  A single character string naming the family.

- dimension:

  The side \\p\\ of the matrix.

- n_free:

  The length \\d\\ of the free vector.

- free_names:

  A character vector of length `n_free`, one label per free value. Fixed
  at construction: every consumer builds parameter tables from these.

- rank:

  The rank of the matrix the family produces. Equal to `dimension` for a
  definite family, and strictly less for a rank-deficient precision,
  where it is the dimension of the space the quadratic form penalises.

- null_basis:

  A `dimension` by `dimension - rank` matrix whose columns are an
  orthonormal basis of the null space, or a matrix with no columns when
  the family is of full rank.

- role:

  One of `"covariance"`, `"precision"` or `"either"`. A label: no method
  reads it and no result depends on it. It exists because the name of a
  family does not say which side of a model it parametrises, and the two
  are different models – the inverse of a compound-symmetry matrix is
  compound symmetry, while the inverse of an AR(1) is tridiagonal and
  not AR(1).

- struct_params:

  A list of whatever the family needs to evaluate itself.

## Value

An object of class `covstruct`. The class is abstract; use one of the
constructors, such as
[`log_cholesky`](https://statmodels7.github.io/covstructs7/reference/log_cholesky.md).

## Details

A structure owns its dimension: `log_cholesky(dimension = 3)` and
`log_cholesky(dimension = 4)` are different objects, with \\d = 6\\ and
\\d = 10\\ free values respectively. The alternative, a dimensionless
recipe applied to whatever arrives, would leave `n_free` and
`free_names` unanswerable before any data exist, and both are needed
then.

The rank and the null space are properties of the family rather than of
a point: scaling a matrix by a positive number and summing positive
semidefinite matrices both leave the null space alone. They are computed
once, at construction, from the components rather than from an assembled
matrix, because the numerical determination of a rank from an assembled
matrix is not scale invariant while the null space is. See
[`struct_null_basis`](https://statmodels7.github.io/covstructs7/reference/struct_null_basis.md).

Only
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md)
is compulsory. Every other generic has a numerical method registered on
this class, so a new structure is a subclass and one method, and a
closed form supplied later replaces the corresponding numerical method
through dispatch.

## See also

[`log_cholesky`](https://statmodels7.github.io/covstructs7/reference/log_cholesky.md),
[`scaled_struct`](https://statmodels7.github.io/covstructs7/reference/scaled_struct.md),
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md),
[`check_covstruct`](https://statmodels7.github.io/covstructs7/reference/check_covstruct.md)

## Examples

``` r
s <- log_cholesky(3)
S7::S7_inherits(s, covstruct)
#> [1] TRUE
c(dimension = s@dimension, n_free = s@n_free, rank = s@rank)
#> dimension    n_free      rank 
#>         3         6         3 
```
