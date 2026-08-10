# Construct a Block Diagonal of Matrix Parameters

\\M(\eta) = \mathrm{diag}(S_1(\eta_1), \ldots, S_B(\eta_B))\\ for matrix
parameters \\S_1, \ldots, S_B\\: a block-diagonal matrix whose blocks
are distinct families, each governed by its own stretch of the free
vector.

## Usage

``` r
block_diag(..., role = NULL)
```

## Arguments

- ...:

  Two or more objects inheriting from
  [`matrix_parameter`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  or a single list of them. Named arguments supply the block labels.

- role:

  One of `"covariance"`, `"precision"` or `"either"`. Defaults to the
  blocks' common role, and to `"either"` when they disagree.

## Value

An object of class
[`BlockDiagParam`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md).

## Details

This is the covariance (or precision) of several independent groups of
coefficients whose structures differ, which is what a model carrying
more than one random-effect term needs. It differs from
[`kron_identity`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
where the blocks are identical and share one free vector; here the free
vectors are concatenated, so the composite has \\\sum_b d_b\\ free
values.

Every quantity of the contract follows from the blocks without
rederivation, and the reason is that the free values of one block do not
enter another: \$\$\partial_k M = \mathrm{diag}(0, \ldots, \partial_k
S_b, \ldots, 0), \qquad \log\lvert M \rvert\_{+} = \sum_b \log\lvert S_b
\rvert\_{+}, \qquad M^{-1} = \mathrm{diag}(S_1^{-1}, \ldots,
S_B^{-1}).\$\$ A derivative whose indices do not all belong to one block
is therefore identically zero, at every order and for the
log-determinant as well as for the value. The rank is the sum of the
blocks' ranks and the null basis is their block diagonal, both read from
the components rather than from an assembled matrix.

The free names are prefixed by the block's label, since two blocks of
the same family would otherwise report the same names and the class
requires them to be unique. Labels come from the names of the arguments
where they are given, and are `b1`, `b2`, ... otherwise.

## See also

[`kron_identity`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)

## Examples

``` r
s <- block_diag(subject = log_cholesky(2), time = ar1(3))
c(dimension = s@dimension, n_free = s@n_free)
#> dimension    n_free 
#>         5         5 
s@free_names
#> [1] "subject_log_L1" "subject_log_L2" "subject_L2.1"   "time_log_scale"
#> [5] "time_z_rho"    
param_logdet(s, c(0, 0, 0, 0, 0))
#> [1] 0
```
