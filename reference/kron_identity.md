# Construct a Block Replication of a Matrix Parameter

Returns an object holding \\M(\eta) = I_m \otimes S(\eta)\\ for an inner
matrix parameter \\S\\ of side \\d\\: a block-diagonal matrix of side
\\md\\ whose \\m\\ diagonal blocks are the **same** matrix, sharing one
free vector.

This is the covariance, or precision, of \\m\\ independent groups whose
within-group structure is common, which is the shape a grouped
random-effect term needs: \\S\\ is the per-group matrix over the
coefficients of one group. Because the blocks share their free values,
`n_free` is the inner parameter's however many groups there are, so a
random slope over a thousand subjects still estimates three numbers.

## Usage

``` r
kron_identity(structure, m)
```

## Arguments

- structure:

  An object inheriting from
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  the per-block parameter. A family that is not a matrix, such as
  [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md),
  throws `'structure' must inherit from 'matrix_parameter'.` The inner
  parameter may itself be rank deficient.

- m:

  The number of blocks, a single integer of at least 1. `0`, a fraction,
  `NA` and a vector all throw
  `'m' must be a single integer of at least 1.` At `m = 1` the result
  equals the inner parameter's value exactly, which makes it safe to
  call in a loop over group counts.

## Value

An object of class
[`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md),
with `dimension` equal to `m * structure@dimension`, `n_free` and
`free_names` the inner parameter's unchanged, `rank` equal to
`m * structure@rank`, `null_basis` the inner one replicated blockwise,
and `param_params` holding `inner` and `m`.

## Every quantity is a linear lift

\$\$\partial_k M = I_m \otimes \partial_k S, \qquad \log\|M\|\_{+} = m
\log\|S\|\_{+}, \qquad M^{-1} = I_m \otimes S^{-1}.\$\$

Nothing is rederived and nothing of side \\md\\ is decomposed: the rank
is \\m\\ times the inner rank, the null basis is the inner one
replicated blockwise, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
loops over the blocks solving the inner parameter \\m\\ times.

The multiplier on the log-determinant is worth expecting. At \\m = 3\\
over a 2 x 2 log-Cholesky covariance,
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
answers `c(6, 6, 0)` where the inner parameter answers `c(2, 2, 0)`: a
step in a shared free value moves every block.

## A deficient inner parameter stays deficient

The rank is exactly \\m\\ times the inner rank, so replicating a
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
on a difference penalty over 2 blocks gives an 8 x 8 matrix of rank 4,
and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
rejects it as the inner one would. The replicated null basis is
annihilated to \\10^{-15}\\.

## The dimension labels

[`kronecker()`](https://rdrr.io/r/base/kronecker.html) drops the inner
parameter's labels, so the value and the four derivative orders are
relabeled `v1`, `v2`, ..., `v(md)` over the composite side. The factor
[`param_factor.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_factor.KronIdentityParam.md)
returns is left bare, as several of the primitive families' are. See
[`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md).

## Notation

\\S\\ is the inner parameter, \\d\\ its side, \\m\\ the number of
blocks, \\I_m\\ the \\m \times m\\ identity and \\\otimes\\ the
Kronecker product. \\\eta\\ is the free vector, shared by every block.

## See also

[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
for blocks that differ and each carry their own free values,
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the other two compositions, and
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
or [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
for the usual inner parameters.

## Examples

``` r
# Three groups sharing one 2 x 2 covariance: 6 x 6, three free values.
inner <- log_cholesky(2)
s <- kron_identity(inner, 3)
c(dimension = s@dimension, n_free = s@n_free)
#> dimension    n_free 
#>         6         3 
eta <- c(0.1, 0, -0.2)
round(param_value(s, eta), 3)
#>        v1     v2     v3     v4     v5     v6
#> v1  1.221 -0.221  0.000  0.000  0.000  0.000
#> v2 -0.221  1.040  0.000  0.000  0.000  0.000
#> v3  0.000  0.000  1.221 -0.221  0.000  0.000
#> v4  0.000  0.000 -0.221  1.040  0.000  0.000
#> v5  0.000  0.000  0.000  0.000  1.221 -0.221
#> v6  0.000  0.000  0.000  0.000 -0.221  1.040

# Every quantity is the inner one lifted, and the log-determinant is
# multiplied by the number of blocks.
c(outer = param_logdet(s, eta), m_times_inner = 3 * param_logdet(inner, eta))
#>         outer m_times_inner 
#>           0.6           0.6 
rbind(outer = param_dlogdet(s, eta), inner = param_dlogdet(inner, eta))
#>       log_L1 log_L2 L2.1
#> outer      6      6    0
#> inner      2      2    0

# The solve is blockwise and exact.
max(abs(param_solve(s, eta) - solve(param_value(s, eta))))
#> [1] 2.775558e-17

# The round trip closes, and m = 1 is the inner parameter itself.
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 6.938894e-17
max(abs(param_value(kron_identity(inner, 1), eta) - param_value(inner, eta)))
#> [1] 0

# A deficient inner parameter stays deficient, block by block.
r <- kron_identity(scaled_matrix(crossprod(diff(diag(4), differences = 2))), 2)
c(dimension = r@dimension, rank = r@rank, null = ncol(r@null_basis))
#> dimension      rank      null 
#>         8         4         4 
max(abs(param_value(r, 0.3) %*% r@null_basis))
#> [1] 1.040818e-15
```
