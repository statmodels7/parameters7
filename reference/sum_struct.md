# Construct a Sum of Fixed Matrices

Returns an object holding \\M(\eta) = \sum_k c_k(\eta_k) P_k\\, for
fixed symmetric positive semidefinite \\P_1, \ldots, P_K\\ and positive
weights carried through a link. The matrices are given once and never
move; what the free vector carries is the \\K\\ weights.

## Usage

``` r
sum_struct(components, link = linkfunctions7::log_link())
```

## Arguments

- components:

  A non-empty list of symmetric positive semidefinite numeric matrices
  of the same side. Each is checked for all three properties, the
  semidefiniteness spectrally at a relative tolerance of \\10^{-8}\\.
  Named entries supply the free-value labels, which must be unique;
  unnamed ones are `w1`, `w2`, ...

- link:

  The positive link carrying each weight onto the free scale,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map onto the positive half line and from the whole
  real line; see
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  for the two conditions.

## Value

An object of class
[`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md),
with `n_free` equal to the number of components, `free_names` the tagged
labels, `null_basis` an orthonormal basis of the components' shared null
space, and `rank` the side less its width.

## What it is for

This is the variance-components covariance \\\sum_k \sigma_k^2
Z_kZ_k^\top\\, and it is also the matrix a penalty with one smoothing
parameter per component assembles. penalties7's `additive_penalty()`
builds the same sum for its own purposes; the difference is that a
penalty is a function of the coefficients while this is a matrix map, so
a distribution can take it as a covariance.

## The value is linear, so most derivative components vanish

\$\$\partial^{m}\_{\eta_k} M = c_k^{(m)}(\eta_k)\\ P_k,\$\$

and a component whose indices name two different weights is **exactly
zero**: at \\K = 2\\ that is 1 of the 3 second-order components, 2 of 4
at third order and 3 of 5 at fourth. What survives is one number times a
fixed matrix.

## The log-determinant is not separable

Its derivatives in the weights come from the cyclic trace expansion

\$\$\frac{\partial^{n}\log\lvert M\rvert} {\partial
c\_{k_1}\cdots\partial c\_{k_n}} = (-1)^{n-1}\sum\_{\sigma}
\operatorname{tr}\bigl(M^{-1}P\_{\sigma(1)}\cdots
M^{-1}P\_{\sigma(n)}\bigr),\$\$

the sum running over the \\(n-1)!\\ cyclic orderings **counted with
multiplicity**, and are then carried onto the free scale by a chain rule
whose Jacobian is diagonal. Counting with multiplicity is load bearing
and the cost of getting it wrong is measured: deduplicating the
orderings that coincide when an index repeats leaves a third-order
component too small by exactly 2 and a fourth-order one by exactly 6.

Against one stencil on the analytic order below, the four orders agree
to \\7 \times 10^{-11}\\, \\1 \times 10^{-11}\\, \\3 \times 10^{-12}\\
and \\3 \times 10^{-13}\\.

## The rank is fixed at construction

The null space of a sum of positive semidefinite matrices is the
**intersection** of theirs, so it does not move with the weights, and it
is read from the components stacked and individually normalized, never
from an assembled matrix. The distinction is measurable, not merely
conceptual. Take \\P_1 = \mathbf{1}\mathbf{1}^\top\\ of side 4, of rank
1, and \\P_2\\ the first-difference penalty, of rank 3: their null
spaces meet only at the origin, so the family has rank 4, and counting
eigenvalues of \\M(\eta)\\ above \\10^{-10}\\ of the largest gives

|  |  |  |  |  |  |
|----|----|----|----|----|----|
| weight ratio | \\10^{0}\\ | \\10^{3.5}\\ | \\10^{6.9}\\ | \\10^{10.4}\\ | \\10^{13.9}\\ |
| eigenvalue count | 4 | 4 | 4 | 1 | 1 |

while the family reports 4 throughout. Weights ten orders of magnitude
apart are an ordinary fitted model. Where there **is** a shared null
space it is held exactly: on the first- and second-difference penalties
over five points, whose null spaces meet in the constants, the residual
\\\lVert M N\rVert / \lVert M\rVert\\ stays at \\3 \times 10^{-16}\\
over fourteen decades of weight ratio.

## What a call costs

The log-determinant's fourth derivatives are the expensive quantity, the
expansion evaluating \\(n-1)!\\ matrix chains per component and one set
partition per repeated index. Seconds per call at side 6, over
repetition loops sized by elapsed time:

|       |               |            |                 |                  |
|-------|---------------|------------|-----------------|------------------|
| \\K\\ | `param_value` | `param_d4` | `param_dlogdet` | `param_d4logdet` |
| 2     | 0.00003       | 0.00017    | 0.00029         | 0.00172          |
| 3     | 0.00002       | 0.00022    | 0.00037         | 0.00344          |
| 5     | 0.00003       | 0.00068    | 0.00039         | 0.01328          |

## Notation

\\K\\ is the number of components, \\P_k\\ the \\k\\-th fixed matrix,
\\c_k = h(\eta_k)\\ its weight, \\p\\ the side, and \\N\\ a basis of the
shared null space.

## See also

[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
for one fixed matrix and one weight,
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
and
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
for the other compositions, and `penalties7::additive_penalty()`, which
assembles the same sum as a penalty.

## Examples

``` r
# Two variance components on three coefficients.
s <- sum_struct(list(between = matrix(1, 3, 3), within = diag(3)))
s@free_names
#> [1] "log_between" "log_within" 
eta <- c(log(0.5), log(2))
param_value(s, eta)
#>     v1  v2  v3
#> v1 2.5 0.5 0.5
#> v2 0.5 2.5 0.5
#> v3 0.5 0.5 2.5

# The value is linear in the weights, so a component naming two of them is
# exactly zero, and a pure one is the weight's derivative times its matrix.
c(mixed = max(abs(param_d2(s, eta)[["log_between:log_within"]])),
  pure = max(abs(param_d1(s, eta)[[1]] - 0.5 * matrix(1, 3, 3))))
#> mixed  pure 
#>     0     0 

# The log-determinant is not separable: the mixed second derivative is as
# large as the pure ones.
param_d2logdet(s, eta)
#> log_between:log_between   log_within:log_within  log_between:log_within 
#>                0.244898                0.244898               -0.244898 

# A single component of rank 1 gives a deficient family, and the
# log-determinant is then the log pseudo-determinant.
d <- sum_struct(list(matrix(1, 3, 3)))
c(rank = d@rank, logdet = param_logdet(d, log(2)), check = log(3 * 2))
#>     rank   logdet    check 
#> 1.000000 1.791759 1.791759 

# The round trip closes.
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 3.330669e-16
```
