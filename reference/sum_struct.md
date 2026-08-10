# Construct a Sum of Fixed Matrices

\\M(\eta) = \sum_k c_k(\eta_k) P_k\\ for fixed symmetric positive
semidefinite \\P_1, \ldots, P_K\\ and positive weights carried through a
link.

## Usage

``` r
sum_struct(
  components,
  link = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- components:

  A list of symmetric matrices of the same side, each positive
  semidefinite. Named entries supply the free-value labels.

- link:

  The positive link carrying each weight onto the free scale. Defaults
  to
  [`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  One of `"covariance"`, `"precision"` or `"either"`.

## Value

An object of class
[`SumStructParam`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md).

## Details

This is the variance-components covariance, \\\sum_k \sigma_k^2
Z_kZ_k'\\, and it is also the matrix a penalty with one smoothing
parameter per component assembles. penalties7's `additive_penalty()`
builds the same sum for its own purposes; the difference is that a
penalty is a function of the coefficients while this is a matrix map, so
a distribution can take it as a covariance.

The value is linear in the weights, so a derivative component is zero
unless every index of the tuple names the same free value, and it then
carries the corresponding derivative of the inverse link times that
component: \$\$\partial^{m}\_{\eta_k} M = c_k^{(m)}(\eta_k)\\ P_k.\$\$
The log-determinant is not separable, and its derivatives in the weights
come from the expansion \$\$\frac{\partial^{n}\log\lvert M\rvert}
{\partial c\_{k_1}\cdots\partial c\_{k_n}} = (-1)^{n-1}\sum\_{\sigma}
\operatorname{tr}\bigl(M^{-1}P\_{\sigma(1)}\cdots
M^{-1}P\_{\sigma(n)}\bigr),\$\$ the sum running over the \\(n-1)!\\
cyclic orderings counted with multiplicity, and are then carried onto
the free scale by a chain rule whose Jacobian is diagonal.

**Rank.** The null space of a sum of positive semidefinite matrices is
the intersection of theirs, so it does not move with the weights and the
rank is fixed at construction. It is read from the components stacked
and individually normalized, never from an assembled matrix: a count of
small eigenvalues of \\M(\eta)\\ falls as the weights spread apart,
which is an ordinary fitted model and not a pathology.

## See also

[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
[`block_diag`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
`additive_penalty`

## Examples

``` r
# two variance components on three coefficients
s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
s@free_names
#> [1] "log_between" "log_within" 
param_value(s, c(log(0.5), log(2)))
#>      [,1] [,2] [,3]
#> [1,]  2.5  0.5  0.5
#> [2,]  0.5  2.5  0.5
#> [3,]  0.5  0.5  2.5
param_logdet(s, c(log(0.5), log(2)))
#> [1] 2.639057
```
