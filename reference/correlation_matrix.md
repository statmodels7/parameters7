# Construct a Correlation Matrix Parameter

A correlation matrix – symmetric positive definite with a unit diagonal
– carried by the spherical parametrization: \\R = LL^\top\\ with each
row of \\L\\ a point on the unit sphere written in angular coordinates,
and each angle carried onto the real line by a bounded link.

## Usage

``` r
correlation_matrix(dimension, role = c("either", "covariance", "precision"))
```

## Arguments

- dimension:

  The side \\p\\ of the matrix.

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Value

An object of class
[`CorrelationParam`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md).

## Details

Row \\i\\ of \\L\\ is built from \\i-1\\ angles \\\theta\_{i1}, \ldots,
\theta\_{i,i-1}\\ in \\(0, \pi)\\: \$\$L\_{ij} = \cos\theta\_{ij}
\prod\_{k \< j} \sin\theta\_{ik} \quad (j \< i), \qquad L\_{ii} =
\prod\_{k \< i} \sin\theta\_{ik}.\$\$ The squared entries of the row
then telescope to one, so \\R = LL^\top\\ has a unit diagonal by
construction rather than by a correction, and it is positive definite
for every value of the angles because \\L\\ is triangular with a
positive diagonal. The construction is that of Rapisarda, Brigo and
Mercurio (2007), and the free values are the angles carried onto
\\\mathbb{R}\\ by
[`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html),
so there is nothing to constrain and no boundary to reach.

Derivatives follow from the same Leibniz rule the log-Cholesky family
uses, since \\R\\ is again a Gram product; what changes is the factor,
whose entries are products of sines and cosines of angles that each
depend on one free value. The rows of \\L\\ are independent, so a
derivative of the *factor* in free values from two different rows
vanishes. The same is not true of \\R\\: its entry \\(i, j)\\ is the
inner product of rows \\i\\ and \\j\\ of \\L\\, so a second derivative
across those two rows need not vanish. What holds instead is that such a
component is supported on exactly the entries \\(i, j)\\ and \\(j, i)\\,
and even there it is zero whenever the two angles differentiated sit
beyond the columns the two rows share.

The free vector runs row by row, and the names `z{i}.{j}` say which row
and which angle, the row.column convention
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
already uses. The ordering is part of the contract: every consumer
builds its parameter tables from these names.

A correlation matrix on its own is what a copula or an LKJ prior is
written against. A covariance is this matrix conjugated by a diagonal of
standard deviations, which is a composition of two parameters rather
than a family of its own.

## References

Rapisarda, F., Brigo, D. and Mercurio, F. (2007). Parameterizing
correlations: a geometric interpretation. *IMA Journal of Management
Mathematics* 18, 55-73.

## See also

[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`compound_symmetry`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`ar1`](https://statmodels7.github.io/parameters7/reference/ar1.md)

## Examples

``` r
s <- correlation_matrix(3)
s@free_names
#> [1] "z2.1" "z3.1" "z3.2"

r <- param_value(s, c(0.4, -0.2, 0.6))
round(r, 4)
#>         v1      v2      v3
#> v1  1.0000 -0.3051  0.1559
#> v2 -0.3051  1.0000 -0.4632
#> v3  0.1559 -0.4632  1.0000
diag(r)
#> v1 v2 v3 
#>  1  1  1 

# the round trip closes exactly
eta <- c(0.4, -0.2, 0.6)
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 2.498002e-16
```
