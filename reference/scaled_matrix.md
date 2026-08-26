# Construct a Scaled Fixed Matrix

Returns an object holding the map \\M(\eta) = h(\eta)\\P\\: a fixed
symmetric positive semidefinite matrix \\P\\, supplied by the caller,
times one positive scale. With `link = NULL` it holds \\P\\ alone and
has no free value.

This is the commonest penalty in semiparametric regression, and it is
the reason the package admits rank-deficient matrices at all. \\P\\ may
be the Gram matrix of a basis derivative, a difference penalty
\\\Delta^\top \Delta\\, a neighborhood matrix, or the identity, which
makes the object a ridge.

## Usage

``` r
scaled_matrix(
  p,
  link = linkfunctions7::log_link(),
  role = c("precision", "covariance", "either"),
  tol = 1e-10
)
```

## Arguments

- p:

  A symmetric positive semidefinite numeric matrix. It must be square,
  free of `NA`, and symmetric to \\10^{-8}\\ relative, and it is
  symmetrized before use. A matrix whose largest eigenvalue is not
  positive throws
  `'p' must be positive semidefinite and not identically zero.`, and one
  whose smallest eigenvalue is below `-tol * max(ev)` throws a message
  quoting both eigenvalues.

- link:

  A linkfunctions7 link carrying the free value onto the positive scale,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map onto the positive half line, so
  `identity_link()` is rejected, and from the whole real line, which
  rules out `sqrt_link()` and its relatives; see
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  for the two conditions. `NULL` means the matrix is fully known:
  `n_free` is then 0, `free_names` is empty, and `param_name` is
  `"fixed"`.

- role:

  A label recording which side of a model the matrix parametrizes,
  defaulting to `"precision"` here because the consumer of a scaled
  fixed matrix is a penalty. No numeric result depends on it.

- tol:

  The relative tolerance below which a singular value of `p` counts as
  zero when its rank is determined, `1e-10` by default. It is passed to
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  and also sets how negative an eigenvalue may be before `p` is
  rejected.

## Value

An object of class
[`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md),
with `n_free` 1 or 0, `free_names` a single link-tagged label or empty,
`rank` and `null_basis` from
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md),
and `param_params` holding `p`, `link` and `logdet_p`.

## Everything is a constant times a function of the scale

Nothing here needs deriving. With the default log link, where \\h(\eta)
= e^{\eta}\\ and \\\lambda = h(\eta)\\,

\$\$\partial\_\eta M = M, \qquad \partial^2\_\eta M = M,\$\$
\$\$\log\|M\|\_+ = r\\\eta + \log\|P\|\_+, \qquad \partial\_\eta
\log\|M\|\_+ = r, \qquad \partial^2\_\eta \log\|M\|\_+ = 0,\$\$

with \\r\\ the rank and \\\log\|P\|\_+\\ the log pseudo-determinant of
\\P\\, computed once at construction and stored. Under another link the
derivatives carry that link's own, through
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md).

## Why the derivative of the log pseudo-determinant matters

It equals the rank, and that is what leaves the scale estimable. Write a
penalty as a negative log prior,

\$\$\tfrac{\lambda}{2}\beta^\top P \beta - \tfrac{r}{2}\log\lambda,\$\$

and the stationary point is \\\lambda = r / (\beta^\top P \beta)\\. Drop
the second term and the derivative keeps one sign, sending the scale to
zero. That second term is the normalizing constant of the prior, which
is why this package keeps it.

## Rank deficiency is admitted, and what it means

A deficient \\P\\ makes the corresponding Gaussian improper, so it is a
legitimate **penalty** without being a legitimate density. Both readings
of a multivariate Gaussian need full rank: a singular covariance is
degenerate on a subspace, and a singular precision does not normalize.
The object still answers
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
with the pseudo-determinant, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
reject it, a consumer of an improper prior needing the quadratic form
and the pseudo-determinant instead of an inverse.

## Notation

\\P\\ is the fixed matrix, \\r\\ its rank, \\\eta\\ the single free
value, \\h = g^{-1}\\ the inverse link and \\\lambda = h(\eta)\\ the
scale. \\\log\|\cdot\|\_+\\ is the log pseudo-determinant, the sum of
the logarithms of the \\r\\ non-zero eigenvalues.

## See also

[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
which is this on the identity from the other direction,
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for a non-negative combination of several fixed matrices,
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for the rank, and
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
for the pseudo-determinant.

## Examples

``` r
# A ridge: the identity, scaled. Full rank, and log|M| = p * eta.
r <- scaled_matrix(diag(4))
c(rank = r@rank, logdet_at_0 = param_logdet(r, 0))
#>        rank logdet_at_0 
#>           4           0 
all.equal(param_logdet(r, 1.3), 4 * 1.3)
#> [1] TRUE

# A second-difference penalty on six coefficients: deficient by two, its
# null space being the constants and the straight lines.
d <- diff(diag(6), differences = 2)
s <- scaled_matrix(crossprod(d))
c(dimension = s@dimension, rank = s@rank)
#> dimension      rank 
#>         6         4 
max(abs(param_value(s, 0.4) %*% s@null_basis))
#> [1] 2.597026e-15

# The derivative of the log pseudo-determinant is the rank, at any scale.
c(at_minus_4 = param_dlogdet(s, -4), at_4 = param_dlogdet(s, 4),
  rank = s@rank)
#> at_minus_4.log_scale       at_4.log_scale                 rank 
#>                    4                    4                    4 

# Which leaves the scale estimable: the stationary point of
# lambda/2 b'Pb - r/2 log(lambda) is r / (b'Pb).
set.seed(2)
b <- rnorm(6)
bPb <- drop(crossprod(b, crossprod(d) %*% b))
obj <- function(lam) lam / 2 * bPb - s@rank / 2 * log(lam)
c(closed_form = s@rank / bPb, numeric = optimize(obj, c(1e-6, 100))$minimum)
#> closed_form     numeric 
#>   0.1250415   0.1250528 

# A deficient family has no inverse and says so.
try(param_solve(s, 0))
#> Error : 'scaled' is rank deficient (4 of 6), so it has no inverse. A consumer of
#>   an improper prior needs the quadratic form and the log
#>   pseudo-determinant, not a pseudo-inverse: assemble the matrix the
#>   model actually inverts and solve that.

# With no link the matrix is fully known and there is nothing to estimate.
f <- scaled_matrix(diag(3), link = NULL)
c(n_free = f@n_free, logdet = param_logdet(f, numeric(0)))
#> n_free logdet 
#>      0      0 
```
