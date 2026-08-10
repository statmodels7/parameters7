# Construct a Covariance as Scales Times a Correlation

\\\Sigma(\eta) = D(\eta_D)\\ R(\eta_R)\\ D(\eta_D)\\ with \\D =
\mathrm{diag}(d_1, \ldots, d_p)\\ carrying the standard deviations
through a positive link, and \\R\\ a correlation matrix parameter.

## Usage

``` r
dr_prod(
  dimension,
  correlation = NULL,
  link = linkfunctions7::log_link(),
  role = c("covariance", "precision", "either")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix, at least 2.

- correlation:

  A
  [`matrix_parameter`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  of side `dimension` producing correlation matrices. Defaults to
  [`correlation_matrix`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)`(dimension)`.

- link:

  The positive link carrying each standard deviation onto the free
  scale. Defaults to
  [`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  One of `"covariance"`, `"precision"` or `"either"`.

## Value

An object of class
[`DrProdParam`](https://statmodels7.github.io/parameters7/reference/DrProdParam.md).

## Details

The separation is what makes the parametrization worth having: the
quantities a reader reads off a fitted covariance are the standard
deviations and the correlations, and here they are the coordinates
rather than a function of them. A log-Cholesky factor produces the same
set of matrices with no coordinate meaning anything on its own.

Every derivative factorizes, because \\\Sigma\_{ij} = d_i d_j R\_{ij}\\
and the two groups of free values are disjoint:
\$\$\partial^{S}\Sigma\_{ij} = \bigl\[\partial^{S_D}(d_i d_j)\bigr\]\\
\bigl\[\partial^{S_R} R\_{ij}\bigr\],\$\$ where \\S_D\\ and \\S_R\\ are
the parts of the multiset \\S\\ falling in each group. The scale factor
is zero unless every index of \\S_D\\ names \\i\\ or \\j\\, so most
components vanish; when \\i = j\\ it is the Leibniz expansion of
\\\partial^m d_i^2\\. Nothing of the correlation family is rederived.

The log-determinant is \\2\sum_j \log d_j + \log\lvert R\rvert\\, which
is separable in the scales and separable from the correlation, so a
log-determinant component mixing two scales, or a scale with a
correlation, is exactly zero.

The correlation block must have full rank. A rank-deficient \\R\\ would
give \\\Sigma\\ the null space \\D^{-1}\ker R\\, which moves with the
free vector, and the class records the rank and the null space as
properties of the family rather than of a point.

## See also

[`correlation_matrix`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`block_diag`](https://statmodels7.github.io/parameters7/reference/block_diag.md)

## Examples

``` r
s <- dr_prod(3)
s@free_names
#> [1] "log_sd1" "log_sd2" "log_sd3" "z2.1"    "z3.1"    "z3.2"   
round(param_value(s, c(0, log(2), 0, 0.4, 0.4, 0.4)), 3)
#>        [,1]   [,2]   [,3]
#> [1,]  1.000 -0.610 -0.305
#> [2,] -0.610  4.000 -0.367
#> [3,] -0.305 -0.367  1.000

# p standard deviations and the p(p-1)/2 angles of the correlation
c(p = 4, n_free = dr_prod(4)@n_free)
#>      p n_free 
#>      4     10 
```
