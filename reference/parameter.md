# Constrained Matrix Parameter

The abstract S7 class of covariance parameters: a map from an
unconstrained vector \\\eta \in \mathbb{R}^{d}\\ to a symmetric matrix
in some constrained set, together with its derivatives and the
quantities a likelihood asks of the matrix.

## Usage

``` r
parameter(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list()
)
```

## Arguments

- param_name:

  A single character string naming the family.

- n_free:

  The length \\d\\ of the free vector.

- free_names:

  A character vector of length `n_free`, one label per free value. Fixed
  at construction: every consumer builds parameter tables from these. A
  label names the coordinate rather than the quantity the coordinate
  produces, and the families here follow one convention for it. Where a
  link carries a constrained quantity onto the free scale, the label
  records that link, so a variance appears as `"log_scale"` and a
  correlation as `"z_rho"`; where the coordinate is already unrestricted
  the label is the plain name of the quantity, as the below-diagonal
  entries `"L2.1"` of a Cholesky factor are. The distinction matters
  outside the family: a consumer flattens the free vector into scalar
  parameters carrying identity links, so a label promising a bounded
  quantity reports a number on a scale that number is not on.

- param_params:

  A list of whatever the family needs to evaluate itself.

## Value

An object of class `parameter`. The class is abstract; use one of the
constructors, such as
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Details

A parameter owns its dimension: `log_cholesky(dimension = 3)` and
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
[`param_null_basis`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md).

Only
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)
is compulsory. Every other generic has a numerical method registered on
this class, so a new parameter is a subclass and one method, and a
closed form supplied later replaces the corresponding numerical method
through dispatch.

## See also

[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md),
[`check_parameter`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)

## Examples

``` r
s <- log_cholesky(3)
S7::S7_inherits(s, parameter)
#> [1] TRUE
c(dimension = s@dimension, n_free = s@n_free, rank = s@rank)
#> dimension    n_free      rank 
#>         3         6         3 
```
