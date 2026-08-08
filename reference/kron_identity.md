# Construct a Block Replication of a Matrix Parameter

\\M(\eta) = I_m \otimes S(\eta)\\ for an inner matrix parameter \\S\\ of
dimension \\d\\: a block-diagonal matrix of dimension \\md\\ whose \\m\\
diagonal blocks are the same matrix, sharing one free vector.

## Usage

``` r
kron_identity(structure, m)
```

## Arguments

- structure:

  An object inheriting from
  [`matrix_parameter`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md):
  the per-block parameter.

- m:

  The number of blocks, a single integer of at least 1.

## Value

An object of class
[`KronIdentityParam`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md).

## Details

This is the covariance (or precision) of \\m\\ independent groups whose
within-group structure is common – the shape a grouped random-effect
term needs, with \\S\\ the per-group matrix over the coefficients of one
group. Every quantity of the contract is a linear lift of the inner
parameter's: \$\$\partial_k M = I_m \otimes \partial_k S, \qquad
\log\lvert M \rvert\_{+} = m \log\lvert S \rvert\_{+}, \qquad M^{-1} =
I_m \otimes S^{-1},\$\$ so nothing is rederived, and the rank is \\m\\
times the inner rank with the null basis replicated blockwise.

The free names are the inner parameter's own: the blocks share their
values by construction, which is what distinguishes this composition
from a general block-diagonal of distinct structures.

## See also

[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)

## Examples

``` r
s <- kron_identity(log_cholesky(2), 3)
c(dimension = s@dimension, n_free = s@n_free)
#> dimension    n_free 
#>         6         3 
param_logdet(s, c(0.1, 0, -0.2))
#> [1] 0.6
```
