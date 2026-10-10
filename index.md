# parameters7

A constrained parameter is a map from an unconstrained vector onto the
set where the parameter lives: the positive definite cone, the simplex,
the row-stochastic matrices.
[parameters7](https://statmodels7.github.io/parameters7/) represents
such a map as an S7 object that can be **evaluated, inverted, and
differentiated to fourth order**. A family whose value is a symmetric
matrix also returns its **log-determinant and its derivatives, its
solves and its factor**. A matrix that is **rank deficient**, as a
penalty precision often is, has a log pseudo-determinant, and its solve
and factor are rejected.

It is the constrained-parameter layer of
[statmodels7](https://statmodels7.github.io), an S7 toolkit for
statistical modeling, alongside
[numericals7](https://statmodels7.github.io/numericals7/),
[linkfunctions7](https://statmodels7.github.io/linkfunctions7/),
[distributions7](https://statmodels7.github.io/distributions7/),
[optimizers7](https://statmodels7.github.io/optimizers7/),
[basis7](https://statmodels7.github.io/basis7/),
[penalties7](https://statmodels7.github.io/penalties7/) and
[modelterms7](https://statmodels7.github.io/modelterms7/).

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/parameters7")
```

Or the whole toolkit at once, which also installs the other seven
packages:

``` r

pak::pak("statmodels7/statmodels7")
```

## A parameter is an object

``` r

s <- log_cholesky(3)
s
#> Parameter: log_cholesky
#> Matrix:    3 x 3, symmetric
#> Rank:      3 of 3
#> 
#> Free values: 6
#>   log_L1, log_L2, log_L3, L2.1, L3.1, L3.2
#> 
#> From the base class: none

eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
round(param_value(s, eta), 4)
#>         v1      v2      v3
#> v1  1.2214  0.5526 -0.4421
#> v2  0.5526  0.9203 -0.0363
#> v3 -0.4421 -0.0363  2.0221
```

The map is a bijection onto the positive definite cone, so the inverse
map recovers the free vector up to rounding:

``` r

max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 6.938894e-17
```

## Parameters and links

A link in
[linkfunctions7](https://statmodels7.github.io/linkfunctions7/) is a
**scalar** bijection with a diagonal Jacobian, and everything built on
it depends on that diagonality. A map from $`\mathbb{R}^d`$ to a matrix
has a full Jacobian, so it is a different class. The two **compose**:
the diagonal of a diagonal matrix is a scalar map, and
[linkfunctions7](https://statmodels7.github.io/linkfunctions7/) links
serve there directly.

``` r

param_value(diagonal_matrix(3, link = linkfunctions7::softplus_link()), c(1, 2, 3))
#>          v1       v2       v3
#> v1 1.313262 0.000000 0.000000
#> v2 0.000000 2.126928 0.000000
#> v3 0.000000 0.000000 3.048587
```

## Beyond covariance matrices

A probability vector lies on the simplex, and each row of a Markov
chain’s transition matrix does too. Both are parameters of the same kind
as a covariance, and their derivatives are written in terms of the value
itself:

``` r

s <- simplex(4)
round(param_value(s, c(0.5, -0.2, 0.8)), 4)
#>     p1     p2     p3     p4 
#> 0.2896 0.1438 0.3909 0.1757

# differentiating sum(pi) = 1 gives zero at every order: each component sums to zero
range(vapply(param_d3(s, c(0.5, -0.2, 0.8)), sum, numeric(1)))
#> [1] -1.387779e-17  3.903128e-18

tm <- transition_matrix(3)
rowSums(param_value(tm, rnorm(6)))
#> s1 s2 s3 
#>  1  1  1
```

The positive definite cone has a second chart, the matrix logarithm
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md).
Its log-determinant is the trace of $`S`$, linear in the free values,
and its inverse is $`\exp(-S)`$. Its derivatives are Fréchet derivatives
of the matrix exponential, computed by the Daleckii–Krein representation
with divided differences that stay accurate under repeated eigenvalues.

``` r

m <- matrix_log(2)
c(logdet = param_logdet(m, c(0.2, -0.3, 0.4)),
  trace = 0.2 - 0.3)
#> logdet  trace 
#>   -0.1   -0.1
```

## Structured covariances

An unstructured covariance over twenty occasions has 210 free values.
The structured families reduce that to two, a scale and a correlation,
with every quantity still in closed form:

``` r

c(unstructured = log_cholesky(20)@n_free,
  compound_symmetry = compound_symmetry(20)@n_free,
  ar1 = ar1(20)@n_free)
#>      unstructured compound_symmetry               ar1 
#>               210                 2                 2

s <- ar1(5)
eta <- c(log(2), atanh(0.6))
round(param_value(s, eta), 3)
#>       v1    v2   v3    v4    v5
#> v1 2.000 1.200 0.72 0.432 0.259
#> v2 1.200 2.000 1.20 0.720 0.432
#> v3 0.720 1.200 2.00 1.200 0.720
#> v4 0.432 0.720 1.20 2.000 1.200
#> v5 0.259 0.432 0.72 1.200 2.000

# the precision of an AR(1) process is tridiagonal, and is returned as such
round(param_solve(s, eta), 3)
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,]  0.781 -0.469  0.000  0.000  0.000
#> [2,] -0.469  1.062 -0.469  0.000  0.000
#> [3,]  0.000 -0.469  1.062 -0.469  0.000
#> [4,]  0.000  0.000 -0.469  1.062 -0.469
#> [5,]  0.000  0.000  0.000 -0.469  0.781
```

An autoregression of any order is
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md).
Its stationary region in the coefficients is not a box (at order two it
is a triangle), so the family is parametrized by the partial
autocorrelations, each in $`(-1, 1)`$ and carried onto the free scale by
the inverse hyperbolic tangent, and the Levinson-Durbin recursion maps
them onto the coefficients. Every free vector is stationary, and the
precision is banded with the order as bandwidth, an autoregression of
order $`q`$ being Markov of order $`q`$:

``` r

s <- autoregressive(7, order = 2)
eta <- c(log(2), atanh(0.7), atanh(-0.3))
round(param_solve(s, eta), 3)
#>        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]
#> [1,]  1.077 -0.980  0.323  0.000  0.000  0.000  0.000
#> [2,] -0.980  1.970 -1.275  0.323  0.000  0.000  0.000
#> [3,]  0.323 -1.275  2.066 -1.275  0.323  0.000  0.000
#> [4,]  0.000  0.323 -1.275  2.066 -1.275  0.323  0.000
#> [5,]  0.000  0.000  0.323 -1.275  2.066 -1.275  0.323
#> [6,]  0.000  0.000  0.000  0.323 -1.275  1.970 -0.980
#> [7,]  0.000  0.000  0.000  0.000  0.323 -0.980  1.077
```

A correlation matrix on its own, as a copula or an LKJ prior needs it,
is
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
whose unit diagonal and positive definiteness follow from the
construction:

``` r

r <- param_value(correlation_matrix(4), c(0.4, -0.2, 0.6, 0.1, -0.5, 0.3))
round(r, 3)
#>        v1     v2     v3     v4
#> v1  1.000 -0.305  0.156 -0.078
#> v2 -0.305  1.000 -0.463  0.380
#> v3  0.156 -0.463  1.000 -0.365
#> v4 -0.078  0.380 -0.365  1.000
```

## Rank-deficient precisions

A spline penalty is singular by construction: its null space contains
the functions that it does not penalize. Read as a prior, the penalty is
improper, so it is a legitimate **penalty** but not a density.

``` r

p <- crossprod(diff(diag(6), differences = 2))
s <- scaled_matrix(p)
s
#> Parameter: scaled
#> Matrix:    6 x 6, symmetric
#> Rank:      4 of 6 (null space of dimension 2)
#> 
#> Free values: 1
#>   log_scale
#> 
#> From the base class: none
```

The null space is spanned by the constants and the straight lines, which
a second-difference penalty does not penalize, and the parameter carries
an orthonormal basis for it:

``` r

lin <- cbind(1, 1:6)
c(penalty_annihilates_them = max(abs(p %*% lin)),
  null_basis_spans_the_same = max(abs(qr.resid(qr(lin), s@null_basis))))
#>  penalty_annihilates_them null_basis_spans_the_same 
#>              0.000000e+00              9.482992e-16
```

The log-determinant becomes the log **pseudo**-determinant, and under
the log link its derivative is the rank at every scale:

``` r

c(at_small_scale = unname(param_dlogdet(s, -6)),
  at_large_scale = unname(param_dlogdet(s, 6)),
  rank = s@rank)
#> at_small_scale at_large_scale           rank 
#>              4              4              4
```

That derivative makes the scale estimable. Writing a penalty as a
negative log prior,
$`\frac{\lambda}{2}\beta^\top P\beta - \frac{r}{2}\log\lambda`$, the
stationary point is $`\lambda = r/(\beta^\top P \beta)`$; without the
term in $`\log\lambda`$ the derivative has one sign and sends
$`\lambda`$ to zero.

``` r

set.seed(1)
beta <- rnorm(6)
q <- drop(t(beta) %*% p %*% beta)

c(closed_form = s@rank / q,
  numerical = optimize(function(l) l / 2 * q - s@rank / 2 * log(l),
                       c(1e-8, 1e8), tol = 1e-12)$minimum)
#> closed_form   numerical 
#>   0.1382685   0.1382685
```

## Rank from the components

The rank of $`\lambda_1 P_1 + \lambda_2 P_2`$ does not depend on the
positive scalars, but a count of eigenvalues above a relative tolerance
does, and a fitted model with smoothing parameters far apart is
ordinary:

``` r

p1 <- kronecker(crossprod(diff(diag(8), differences = 2)), diag(4))
p2 <- kronecker(diag(8), crossprod(diff(diag(4), differences = 2)))

count_ev <- function(m) {
  e <- svd(m, nu = 0, nv = 0)$d
  sum(e > 1e-8 * max(e))
}

c(balanced = count_ev(p1 + p2),
  ratio_1e12 = count_ev(1e-6 * p1 + 1e6 * p2),
  from_components = param_null_basis(list(p1, p2))$rank)
#>        balanced      ratio_1e12 from_components 
#>              28              16              28
```

[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
therefore stacks the individually normalized components (the null space
of a sum of positive semidefinite matrices is the intersection of their
null spaces), and membership is tested against that basis.

## Consumers of a parameter

A multivariate distribution in
[distributions7](https://statmodels7.github.io/distributions7/) carries
its covariance as a parameter, on either side: a parameter imposed on
$`\Sigma`$ and the same parameter imposed on $`\Omega = \Sigma^{-1}`$
are different models whenever the family is not closed under inversion.
The free values become parameters of the distribution, so the
derivatives, the validator and the fitting routine of that package apply
with no special case.

``` r

# in distributions7; not run here, since this package does not depend on its
# own consumer
d <- mvgaussian1_distrib(2, parameters7::log_cholesky(2))
d@params
#> [1] "mu1" "mu2" "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"

fit <- fit_distrib(d, y)
mv_sigma(d, coef(fit))
```

The names above are the parameter’s own `free_names` with a prefix
naming the matrix that they build. They are fixed at construction
because every consumer builds its parameter tables from them.

A consumer reports interpretable quantities, not coordinates.
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
returns the quantities that a family declares (a marginal variance, a
correlation, the coefficients of an autoregression) together with the
Jacobian of the map from the free vector and the scale on which each
interval is built, so that the consumer can report them by the delta
method.

``` r

p <- param_readable(autoregressive(8, order = 2), c(log(4), atanh(0.75), atanh(-0.35)))
round(p$value, 4)
#>   scale   pacf1   pacf2    phi1    phi2 
#>  4.0000  0.7500 -0.3500  1.0125 -0.3500
p$transform
#>      scale      pacf1      pacf2       phi1       phi2 
#>      "log"    "atanh"    "atanh" "identity" "identity"
```

## A user-defined parameter

The derivatives, the log-determinant, the solve and the factor have
methods on the base classes, so a new parameter is a subclass with one
method,
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).
The inverse map
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
has no default method.

``` r

ExpDecay <- S7::new_class("ExpDecay", parent = matrix_parameter)

S7::method(param_value, ExpDecay) <- function(s, eta, ...) {
  d <- abs(outer(s@param_params$times, s@param_params$times, "-"))
  m <- exp(eta[1]) * exp(-d / exp(eta[2]))
  dimnames(m) <- rep(list(paste0("v", seq_len(nrow(m)))), 2)
  m
}

decay <- ExpDecay(
  param_name = "exp_decay", dimension = 4L, n_free = 2L,
  free_names = c("log_scale", "log_range"), rank = 4L,
  null_basis = matrix(numeric(0), 4, 0),
  param_params = list(times = c(0, 0.5, 1.7, 3))
)

# not implemented by the family, computed by the base methods
round(param_d1(decay, c(0.2, 0.6))[["log_range"]], 4)
#>        v1     v2     v3     v4
#> v1 0.0000 0.2547 0.4483 0.3876
#> v2 0.2547 0.0000 0.4163 0.4250
#> v3 0.4483 0.4163 0.0000 0.4269
#> v4 0.3876 0.4250 0.4269 0.0000
param_logdet(decay, c(0.2, 0.6))
#> [1] -0.6482252
```

The example generalizes an AR(1) correlation to unequally spaced times,
which
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
does not cover. Its free names follow the convention of the package: a
name identifies the coordinate and not the quantity that the coordinate
produces. Where a link carries a constrained quantity onto the free
scale the name records that link, so a scale appears as `log_scale` and
the correlation of
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) as
`z_rho`; where the coordinate is already unrestricted, as the
below-diagonal entries of a Cholesky factor are, the name is the plain
`L2.1`. A consumer flattens the free vector into scalars carrying
identity links, so a name implying a bounded quantity would show a
number on a scale to which it does not belong.

## Validating a parameter

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs nine checks on a matrix family, fewer on a family that is not a
matrix, each against a route that the implementation does not take. A
quantity that comes from the base class is reported as **not checked**:
comparing a finite difference with a finite difference is the same
arithmetic twice.

``` r

invisible(check_parameter(log_cholesky(3)))
#> Parameter: log_cholesky   (3 x 3, rank 3, 6 free)
#>   [OK         ] membership           0.00e+00
#>   [OK         ] round trip           7.77e-16
#>   [OK         ] first derivatives    5.16e-11
#>   [OK         ] second derivatives   4.45e-11
#>   [OK         ] log-determinant      6.13e-15
#>   [OK         ] logdet gradient      1.62e-13
#>   [OK         ] logdet hessian       0.00e+00
#>   [OK         ] solve and factor     5.80e-15
#>   [OK         ] shapes and names  
#>   9 passed, 0 failed, 0 not checked
invisible(check_parameter(decay))
#> Parameter: exp_decay   (4 x 4, rank 4, 2 free)
#>   [OK         ] membership           0.00e+00
#>   [NOT CHECKED] round trip        
#>   [NOT CHECKED] first derivatives 
#>   [NOT CHECKED] second derivatives
#>   [NOT CHECKED] log-determinant   
#>   [NOT CHECKED] logdet gradient   
#>   [NOT CHECKED] logdet hessian    
#>   [OK         ] solve and factor     6.41e-16
#>   [OK         ] shapes and names  
#>   3 passed, 0 failed, 6 not checked
```

## Contents

|  |  |
|----|----|
| unstructured matrices | [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md), [`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md), [`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md) |
| structured matrices | [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md), [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md), [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md), [`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md), [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md) |
| diagonal and fixed | [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md), [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md), [`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md) |
| compositions | [`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md), [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md), [`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md), [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md), [`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md) |
| not symmetric matrices | [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md), [`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md) |
| maps | [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md), [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md) |
| derivatives | [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md) … [`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md), exact to fourth order |
| likelihood pieces | [`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md), [`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md) … [`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md), [`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md), [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md), [`param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md), [`param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md) |
| interpretable quantities | [`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md) |
| tools | [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md), [`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md), [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md), [`print()`](https://rdrr.io/r/base/print.html) |

## Related

- [numericals7](https://statmodels7.github.io/numericals7/) — the
  numerical layer of the toolkit: stencils, quadrature, special
  functions
- [linkfunctions7](https://statmodels7.github.io/linkfunctions7/) — link
  functions with exact derivatives to fifth order
- [distributions7](https://statmodels7.github.io/distributions7/) —
  distributions carrying exact derivatives of the log-likelihood
- [optimizers7](https://statmodels7.github.io/optimizers7/) —
  optimization algorithms and stopping rules as objects
- [basis7](https://statmodels7.github.io/basis7/) — basis expansions
  with exact derivatives, integrals and Gram matrices
- [penalties7](https://statmodels7.github.io/penalties7/) — penalties as
  objects
- [modelterms7](https://statmodels7.github.io/modelterms7/) — model
  terms as objects
- [the book](https://statmodels7.github.io/book/) — the mathematics
  behind the toolkit
