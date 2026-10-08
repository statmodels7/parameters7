# Construct the Inverse of a Matrix Parameter

Returns the family \\N(\eta) = S(\eta)^{-1}\\ for a matrix family \\S\\:
the same free vector, the inverted value. It is the composition that
says "the other side" to a consumer that fixes one, as
`penalties7::structured_penalty()` fixes the precision.

## Usage

``` r
inverse_of(structure)
```

## Arguments

- structure:

  An object inheriting from
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  of full rank.

## Value

An object of class
[`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md),
with `dimension`, `rank`, `n_free` and `free_names` the inner family's,
`null_basis` a `dimension` by 0 matrix, and `param_name`
`"inverse_of(inner)"`.

## What it is for

A family and its inverse are the same set of matrices only when the
family is closed under inversion, which
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
and
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
are and
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
are not. For a family that is not closed the two sides are different
models, and this wrapper is how the one the family does not name is
written: `inverse_of(correlation_matrix(3))` is the matrix whose inverse
is a correlation matrix, which no chart of the package produces
directly.

Where the inverse has a structure worth exploiting the package writes it
out instead:
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
and
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
give the same values as `inverse_of(ar1())` and
`inverse_of(autoregressive())` and reach them in \\O(p)\\ entries rather
than through a factorization.

## The derivatives

Differentiating \\N = S^{-1}\\ repeatedly gives a sum over the
**ordered** set partitions of the differentiated positions,

\$\$\partial_I N = \sum\_{(B_1, \ldots, B_q)} (-1)^q\\
N\\(\partial\_{B_1} S)\\N \cdots N\\(\partial\_{B_q} S)\\N,\$\$

the blocks being ordered because matrices do not commute. At first order
it is \\-N (\partial_k S) N\\ and at second \\N(\partial_k
S\\N\\\partial_l S + \partial_l S\\N\\\partial_k S)N -
N(\partial\_{kl}S)N\\, which are the two expressions written by hand
elsewhere in the toolkit. The number of terms is the Fubini number of
the order: 1, 3, 13 and 75.

Every factor is the inner family's own closed form, so nothing is
differenced: the arithmetic is exact wherever the inner family's is.

## The log-determinant

\\\log\lvert N\rvert = -\log\lvert S\rvert\\, so the value and all four
derivative orders are the inner family's negated, and no determinant of
the inverted matrix is taken.

## Rank

The inner family must be full rank. A singular matrix has no inverse, so
there is no family to build, and the constructor rejects it naming the
rank and the dimension.

## Notation

\\S\\ is the inner family's map, \\N = S^{-1}\\ the one this builds,
\\I\\ a multiset of differentiated positions and \\(B_1,\ldots,B_q)\\ an
ordered partition of it into non-empty blocks.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
and
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
for the two written-out inverses,
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the other compositions, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md),
which a family with a closed inverse implements and which this reads.

## Examples

``` r
# The inverse of an AR(1) correlation is tridiagonal, so the family whose
# value is that inverse is a different model from ar1() itself.
s <- inverse_of(ar1(5))
eta <- c(log(2), atanh(0.6))
round(param_value(s, eta), 4)
#>         v1      v2      v3      v4      v5
#> v1  0.7812 -0.4688  0.0000  0.0000  0.0000
#> v2 -0.4688  1.0625 -0.4688  0.0000  0.0000
#> v3  0.0000 -0.4688  1.0625 -0.4688  0.0000
#> v4  0.0000  0.0000 -0.4688  1.0625 -0.4688
#> v5  0.0000  0.0000  0.0000 -0.4688  0.7812

# It is the inner family's inverse, exactly.
max(abs(param_value(s, eta) - solve(param_value(ar1(5), eta))))
#> [1] 2.220446e-16

# The log-determinant is the inner one negated.
c(inverse = param_logdet(s, eta), inner = param_logdet(ar1(5), eta))
#>   inverse     inner 
#> -1.680587  1.680587 

# The first derivative is -N (d S) N.
N <- param_value(s, eta)
A <- param_d1(ar1(5), eta)
max(abs(param_d1(s, eta)[[1]] - (-N %*% A[[1]] %*% N)))
#> [1] 0

# The free vector is the inner family's, so the round trip closes on it.
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 3.330669e-16
```
