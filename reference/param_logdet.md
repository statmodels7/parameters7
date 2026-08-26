# Log-Determinant of a Parameter's Matrix

Returns \\\log\|M\|\\ for a full-rank family, and the log
pseudo-determinant, the sum of the logs of the non-zero eigenvalues, for
a rank-deficient one. This is the quantity a Gaussian likelihood needs
beside the quadratic form, and every matrix family answers it, in closed
form where one exists.

One generic covers both cases because a consumer asks the same question
of either: what normalizing constant does this matrix contribute. Which
answer is the right one follows from the object's declared `rank`, so
the caller does not have to branch.

## Usage

``` r
param_logdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).
  A family whose value is not a symmetric matrix has no method, and the
  call fails at dispatch with
  `Can't find method for param_logdet(<SimplexParam>)`.

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A single number. `-Inf` is not returned: a full-rank family is positive
definite at every finite `eta`, and a deficient one drops its null
directions instead of taking `log(0)`.

## The sign belongs to the caller

A Gaussian log-density written in the covariance carries
\\-\tfrac{1}{2}\log\|\Sigma\|\\ and one written in the precision carries
\\+\tfrac{1}{2}\log\|\Omega\|\\. The parameter reports what the
log-determinant is; the likelihood decides where it goes and with what
sign. The object's `role` records which side it was built for, and no
method reads it.

## Computed from the parametrization, not from the matrix

A closed form is usually far cheaper than a decomposition, and often
simply linear. In the log-Cholesky parametrization \\\log\|M\| =
2\sum\_{i=1}^{p} \eta_i\\, twice the sum of the first \\p\\ free values,
so no factorization happens at all. The
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
family answers \\p\\\eta_1 + (p-1)\log(1 - \rho^2)\\, again in constant
work whatever \\p\\ is. The base method on
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
which a family written elsewhere inherits, takes an eigendecomposition
and sums the logs of the eigenvalues above a relative tolerance, at
\\O(p^3)\\.

## Rank deficiency

A deficient family returns the log pseudo-determinant, and its
`null_basis` says which directions were left out. That is the quantity
an improper prior contributes to a marginal likelihood: the penalized
normal equations invert \\X^\top X + \lambda P\\, which is non-singular
even when \\P\\ is not, so the deficiency never has to be inverted.
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
reject it for that reason.

## Notation

\\M\\ is the matrix
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns, \\\Sigma\\ read as a covariance and \\\Omega\\ read as a
precision. \\\eta\\ is the free vector and \\p\\ the side of the matrix.

## See also

[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for its derivatives,
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the quadratic form that goes with it, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for the directions a deficient family omits.

## Examples

``` r
# It agrees with the eigenvalues, and is cheaper: no decomposition is taken.
s <- log_cholesky(3)
eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
M <- param_value(s, eta)
c(param_logdet(s, eta), determinant(M, logarithm = TRUE)$modulus)
#> [1] 0.4 0.4

# In this parametrization it is linear: twice the sum of the first p entries.
all.equal(param_logdet(s, eta), 2 * sum(eta[1:3]))
#> [1] TRUE

# An AR(1) covariance, where it is not linear.
a <- ar1(4)
a@free_names
#> [1] "log_scale" "z_rho"    
eta_a <- c(0.3, 0.8)
rho <- tanh(eta_a[2])
all.equal(param_logdet(a, eta_a), 4 * eta_a[1] + 3 * log(1 - rho^2))
#> [1] TRUE

# A rank-deficient precision: the pseudo-determinant over the four non-zero
# eigenvalues of a second-difference penalty on six coefficients.
P <- crossprod(diff(diag(6), differences = 2))
r <- scaled_matrix(P)
ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
c(rank = r@rank, logdet = param_logdet(r, 0), from_ev = sum(log(ev[1:4])))
#>    rank  logdet from_ev 
#> 4.00000 4.65396 4.65396 

# The scale enters it rank times over, so a step of 1 in the free value
# moves the log pseudo-determinant by the rank.
param_logdet(r, 1) - param_logdet(r, 0)
#> [1] 4
```
