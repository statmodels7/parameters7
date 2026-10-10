# Construct a Covariance as Scales Times a Correlation

Returns an object holding \\\Sigma(\eta) = D(\eta_D)\\ R(\eta_R)\\
D(\eta_D)\\, with \\D = \mathrm{diag}(d_1, \ldots, d_p)\\ carrying the
standard deviations through a positive link and \\R\\ a correlation
matrix parameter.

## Usage

``` r
dr_prod(dimension, correlation = NULL, link = linkfunctions7::log_link())
```

## Arguments

- dimension:

  The side \\p\\ of the matrix, at least 2: a 1 by 1 correlation matrix
  has no free value, and the constructor signals an error.

- correlation:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  of side `dimension` **producing correlation matrices**, that is with a
  unit diagonal at every free vector. Defaults to
  [`correlation_matrix(dimension)`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md);
  [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
  and
  [`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
  of
  [`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
  blocks have the property as well. A block carrying a scale of its own,
  such as
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
  makes the composite unidentified, the same matrix arising from a whole
  ray of free vectors; the constructor reads the diagonal at two probe
  free vectors and rejects such a block, which
  [`check_unit_diagonal()`](https://statmodels7.github.io/parameters7/reference/check_unit_diagonal.md)
  describes.

- link:

  The positive link carrying each standard deviation onto the free
  scale,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map onto the positive half line and from the whole
  real line; see
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  for the two conditions.

## Value

An object of class
[`DrProdParam()`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md),
with `n_free` equal to \\p\\ plus the correlation's, `free_names` the
tagged `log_sd1` ... `log_sdp` followed by the correlation's own, `rank`
equal to `dimension`, and an empty `null_basis`.

## The separation of scales and correlations

The quantities that a reader takes off a fitted covariance are the
standard deviations and the correlations. Here the standard deviations
are the inverse link of the first \\p\\ coordinates, and the
correlations are the value of the correlation block at the remaining
ones. A
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
factor produces the same set of matrices with coordinates that have no
separate meaning. The example below reads the standard deviations off
the diagonal and recovers the correlation block's value by dividing them
out.

## Derivatives

Because \\\Sigma\_{ij} = d_i d_j R\_{ij}\\ and the two groups of free
values are disjoint,

\$\$\partial^{S}\Sigma\_{ij} = \bigl\[\partial^{S_D}(d_i d_j)\bigr\]\\
\bigl\[\partial^{S_R} R\_{ij}\bigr\],\$\$

where \\S_D\\ and \\S_R\\ are the parts of the multiset \\S\\ falling in
each group. The correlation family's own components are fetched and
multiplied entrywise, with no new derivation.

Both factors are **sparse**, and where their supports miss each other
the component is exactly zero. The scale factor is supported on the rows
and columns that the indices of \\S_D\\ name, and vanishes altogether
once \\S_D\\ names three distinct scales; a correlation derivative in an
angle of row \\i\\ is supported on the entries of row and column \\i\\
whose inner product involves that angle, which for `z3.2` are \\(3,2)\\
and \\(2,3)\\. So `log_sd1:z3.2` is zero, the first factor living in row
and column 1 and the second on entries \\(3,2)\\ and \\(2,3)\\.

## The log-determinant

\$\$\log\lvert\Sigma\rvert = 2\sum_j \log d_j + \log\lvert R\rvert,\$\$

separable in the scales and separable from the correlation, so a
component mixing two scales, or a scale with a correlation, is exactly
zero. Under the default log link the first derivative in a scale is 2 at
every point, the scales entering on both sides; for another link it is
\\2h'(\eta)/h(\eta)\\.

## The correlation block must have full rank

A rank-deficient \\R\\ would give \\\Sigma\\ the null space \\D^{-1}\ker
R\\, which **moves with the free vector**, while the class records the
rank and the null space as properties of the family. The constructor
rejects one instead of recording a rank that is not stable.

## Notation

\\p\\ is the dimension, \\d_j\\ the \\j\\-th standard deviation, \\D =
\mathrm{diag}(d)\\, \\R\\ the correlation matrix, and \\\eta_D\\ and
\\\eta_R\\ the two stretches of the free vector.

## See also

[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
for the default block,
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for the same set of matrices in coordinates without a separate meaning,
and
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the other compositions.

## Examples

``` r
s <- dr_prod(3)
s@free_names
#> [1] "log_sd1" "log_sd2" "log_sd3" "z2.1"    "z3.1"    "z3.2"   

eta <- c(log(1), log(2), log(0.5), 1.0, 1.2, 0.9)
M <- param_value(s, eta)

# The standard deviations are the inverse link of the first coordinates,
# read back off the diagonal.
rbind(from_matrix = sqrt(diag(M)), from_eta = exp(eta[1:3]))
#>             v1 v2  v3
#> from_matrix  1  2 0.5
#> from_eta     1  2 0.5

# And dividing them out leaves the correlation block's own value.
R <- M / outer(sqrt(diag(M)), sqrt(diag(M)))
max(abs(R - param_value(correlation_matrix(3), eta[4:6])))
#> [1] 0

# A second derivative is exactly zero where the two factors' supports miss
# each other: the scale factor lives in row and column 1, the correlation
# derivative on the (3, 2) entry.
max(abs(param_d2(s, eta)[["log_sd1:z3.2"]]))
#> [1] 0

# The log-determinant separates, so the first derivative in a scale is 2.
param_dlogdet(s, eta)[1:3]
#> log_sd1 log_sd2 log_sd3 
#>       2       2       2 

# The round trip closes.
max(abs(param_free(s, M) - eta))
#> [1] 6.661338e-16
```
