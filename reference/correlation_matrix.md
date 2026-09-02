# Construct a Correlation Matrix Parameter

Returns an object holding the spherical parametrization of a correlation
matrix: symmetric, positive definite, and with a unit diagonal that
holds exactly, with no correction applied afterwards. Each row of the
Cholesky factor is a point on the unit sphere written in angular
coordinates, and each angle is carried onto the whole real line by a
bounded link, so any free vector in \\\mathbb{R}^{p(p-1)/2}\\ gives a
valid correlation matrix.

Use it when the correlations are what the model is about and a scale
belongs elsewhere: a copula, an LKJ prior, or the correlation half of a
covariance built with
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md).

## Usage

``` r
correlation_matrix(dimension)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix. A single positive whole number, finite
  and at least 1; \\p = 1\\ gives a constant with no free values, as
  **Details** describes. Anything else throws
  `'dimension' must be a single positive integer.`

## Value

An object of class
[`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md),
with `n_free` equal to \\p(p-1)/2\\, `free_names` `z2.1`, `z3.1`,
`z3.2`, ... row by row, `rank` equal to `dimension`, an empty
`null_basis`, `param_name` `"correlation"`, and `param_params` holding
`row`, `col` and `link`.

## The construction

Row \\i\\ of \\L\\ is built from \\i-1\\ angles \\\theta\_{i1}, \dots,
\theta\_{i,i-1}\\ in \\(0, \pi)\\:

\$\$L\_{ij} = \cos\theta\_{ij} \prod\_{k \< j} \sin\theta\_{ik} \quad (j
\< i), \qquad L\_{ii} = \prod\_{k \< i} \sin\theta\_{ik}.\$\$

The squared entries of a row telescope to one, so every row of \\L\\ is
a unit vector and \\R = LL^\top\\ has a unit diagonal by construction.
It is positive definite at every value of the angles, \\L\\ being
triangular with a positive diagonal. The angles reach the free scale
through `linkfunctions7::bounded_link(lwr = 0, upr = pi)`, so there is
nothing to constrain and no boundary to run into.

## How the derivatives behave, and a claim that does not hold

Derivatives come from the same Leibniz rule the log-Cholesky family
uses, \\R\\ being a Gram product again; what changes is the factor,
whose entries are products of sines and cosines of angles that each
depend on one free value.

The rows of \\L\\ are independent, so a derivative of the **factor** in
free values from two different rows is zero. That is **not** true of
\\R\\: its entry \\(i, j)\\ is the inner product of rows \\i\\ and \\j\\
of \\L\\, so a second derivative across those two rows need not vanish.
What does hold is that such a component is supported on exactly the
entries \\(i, j)\\ and \\(j, i)\\, and is zero even there when the two
angles differentiated sit beyond the columns the two rows share. A
structural claim about a product is not a claim about its factors.

## The log-determinant is separable

The factor is triangular, so \\\|R\| = \prod_i L\_{ii}^2\\ and

\$\$\log\|R\| = 2 \sum\_{i,k} \log \sin\theta\_{ik},\$\$

one term per free value. Every mixed derivative of the log-determinant
is therefore exactly zero at every order, and each pure one is a
logarithm composed with that angle's sine. No factorization and no
determinant is computed.

## Reading the free vector

It runs row by row, and the names `z{i}.{j}` say which row and which
angle, the row.column convention
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
uses. The `z` records the link, so a free value of 0.4 is not an angle
of 0.4; the angle is `bounded_link(0, pi)`'s inverse of it. The ordering
is part of the interface.

## Where double precision gives out

Measured at \\p = 3\\ with every free value equal: at 5 the smallest
eigenvalue is \\3 \times 10^{-8}\\ and at 10 it is \\-2 \times
10^{-16}\\, the matrix having reached its boundary in double precision
with a correlation of \\-1\\ to every printed digit. The parametrization
is exact; what gives out is the representation of a correlation one ulp
from the edge. A fit that walks a free value past about \\\pm 8\\ is
reporting a boundary, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
sweeps to \\\pm 2\\ for that reason.

## p = 1 is a constant

A \\1 \times 1\\ correlation matrix has no angles, so `n_free` is 0 and
`free_names` is empty.
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns the \\1 \times 1\\ identity at the only free vector there is,
`numeric(0)`,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
return empty lists, and
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
is 0.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
passes, reporting its two log-determinant derivative rows as
`NOT CHECKED`, there being no free value to differentiate in. The family
therefore degenerates to a constant rather than refusing, so it stays
composable inside
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md);
it carries no information of its own.

## Notation

\\\eta\\ is the free vector, of length \\d = p(p-1)/2\\, and \\p\\ the
side of the matrix. \\\theta\_{ij} \in (0, \pi)\\ is the \\j\\-th angle
of row \\i\\, \\L\\ the lower triangular factor and \\R = LL^\top\\ the
correlation matrix. Note that \\\theta\\ is an angle here and a
distribution parameter elsewhere in the toolkit.

## References

Rapisarda, F., Brigo, D. and Mercurio, F. (2007). Parameterizing
correlations: a geometric interpretation. *IMA Journal of Management
Mathematics* **18**, 55-73.

## See also

[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
to combine this with a diagonal of standard deviations into a
covariance,
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for an unstructured matrix,
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
and
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
for structured correlations with two free values, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
for the map.

## Examples

``` r
# Three angles for a 3 x 3 correlation matrix, named by row and column.
s <- correlation_matrix(3)
s@free_names
#> [1] "z2.1" "z3.1" "z3.2"

eta <- c(0.4, -0.2, 0.6)
r <- param_value(s, eta)
round(r, 4)
#>         v1      v2      v3
#> v1  1.0000 -0.3051  0.1559
#> v2 -0.3051  1.0000 -0.4632
#> v3  0.1559 -0.4632  1.0000

# The unit diagonal is exact, and the matrix is positive definite.
diag(r)
#> v1 v2 v3 
#>  1  1  1 
eigen(r, only.values = TRUE)$values > 0
#> [1] TRUE TRUE TRUE

# The round trip closes exactly.
max(abs(param_free(s, r) - eta))
#> [1] 2.498002e-16

# Absurd free values stay in the set, which is the point.
m <- param_value(s, c(-30, 40, -20))
c(min_diag = min(diag(m)), max_abs_corr = max(abs(m[upper.tri(m)])))
#>     min_diag max_abs_corr 
#>            1            1 

# The log-determinant is twice the sum of the logs of the sines, so it needs
# no factorization, and it agrees with the eigenvalues.
theta <- linkfunctions7::linkinv(s@param_params$link, eta)
c(closed_form = param_logdet(s, eta), from_sines = 2 * sum(log(sin(theta))),
  from_eigen = sum(log(eigen(r, only.values = TRUE)$values)))
#> closed_form  from_sines  from_eigen 
#>  -0.3394489  -0.3394489  -0.3394489 

# And it is separable, so every mixed second derivative is exactly zero.
round(param_d2logdet(s, eta), 10)
#> z2.1:z2.1 z3.1:z3.1 z3.2:z3.2 z2.1:z3.1 z2.1:z3.2 z3.1:z3.2 
#> -1.160942 -1.214977 -1.077535  0.000000  0.000000  0.000000 
```
