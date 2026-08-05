# Construct a Scaled Fixed Matrix

\\M(\eta) = h(\eta)\\P\\ for a fixed symmetric positive semidefinite
\\P\\ and a positive scale \\h(\eta)\\, or the fixed \\P\\ itself when
no link is given.

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

  A symmetric positive semidefinite matrix.

- link:

  A linkfunctions7 link mapping the free scale to the positive scale, or
  `NULL` for a fully known matrix with no free value. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).
  Defaults to `"precision"`, since the consumer of a scaled fixed matrix
  is a penalty.

- tol:

  The relative tolerance below which a singular value of `p` counts as
  zero when its rank is determined.

## Value

An object of class
[`ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md).

## Details

This is the commonest penalty in semiparametric regression, and the
reason the package admits rank-deficient matrices at all. \\P\\ may be a
Gram matrix of a basis derivative, a difference penalty \\\Delta^\top
\Delta\\, a neighborhood matrix, or the identity, which makes the
parameter a ridge.

Everything is a constant times a function of the scale, so nothing needs
deriving. With the default log link, where \\h(\eta) = e^{\eta}\\ and
\\\lambda = h(\eta)\\, \$\$\partial\_\eta M = M, \qquad \partial^2\_\eta
M = M,\$\$ \$\$\log\|M\|\_+ = r\\\eta + \log\|P\|\_+, \qquad
\partial\_\eta \log\|M\|\_+ = r, \qquad \partial^2\_\eta \log\|M\|\_+ =
0,\$\$ with \\r\\ the rank and \\\log\|P\|\_+\\ computed once at
construction.

The derivative of the log pseudo-determinant being the rank is what
makes the scale estimable. Writing a penalty as a negative log prior,
\\\frac{\lambda}{2}\beta^\top P \beta - \frac{r}{2}\log\lambda\\, the
stationary point is \\\lambda = r / (\beta^\top P \beta)\\, whereas
dropping the second term leaves a derivative of one sign and sends the
scale to zero. That second term is exactly the normalizing constant of
the prior, which is the reason this package keeps it.

A rank-deficient \\P\\ makes the corresponding gaussian improper, which
is what makes it a legitimate penalty and not a legitimate density. Both
readings of a multivariate gaussian require full rank: a singular
covariance is degenerate on a subspace, and a singular precision does
not normalize.

## See also

[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`param_null_basis`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)

## Examples

``` r
# a ridge: the identity, scaled
r <- scaled_matrix(diag(4))
c(rank = r@rank, logdet = param_logdet(r, 0))
#>   rank logdet 
#>      4      0 

# a second-difference penalty: rank deficient, and its null space is the
# constants and the straight lines, which the penalty leaves alone
d <- diff(diag(6), differences = 2)
s <- scaled_matrix(crossprod(d))
c(dimension = s@dimension, rank = s@rank)
#> dimension      rank 
#>         6         4 

# the derivative of the log pseudo-determinant is the rank, at any scale
c(param_dlogdet(s, -4), param_dlogdet(s, 4))
#> log_scale log_scale 
#>         4         4 

# and a fully known matrix has no free value at all
scaled_matrix(diag(3), link = NULL)@n_free
#> [1] 0
```
