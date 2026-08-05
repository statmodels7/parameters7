# Construct an Autoregressive Parameter

The covariance of \\p\\ consecutive observations of a stationary
autoregression of order \\q\\, \$\$y_t = \phi_1 y\_{t-1} + \cdots +
\phi_q y\_{t-q} + \varepsilon_t,\$\$ parametrized by its marginal
variance and its \\q\\ partial autocorrelations, so \\q + 1\\ free
values whatever the dimension.

## Usage

``` r
autoregressive(
  dimension,
  order,
  link_scale = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix: the number of consecutive observations.
  Must exceed the order, since a shorter stretch does not identify the
  last partial autocorrelation.

- order:

  The order \\q\\ of the autoregression, at least 1.

- link_scale:

  A linkfunctions7 link onto the positive scale, carrying the marginal
  variance. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Value

An object of class
[`AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md).

## Details

The parametrization is forced by the shape of the stationary region. The
coefficients \\\phi\\ are stationary exactly when the roots of \\1 -
\phi_1 z - \cdots - \phi_q z^{q}\\ lie outside the unit circle, and that
set is not a box: already at \\q = 2\\ it is a triangle, so no
collection of scalar links onto intervals can cover it. The partial
autocorrelations do not have this problem. Each lies in \\(-1, 1)\\
independently of the others, and the Levinson-Durbin recursion carries
them onto the stationary coefficients bijectively, which is the
transformation of Barndorff-Nielsen and Schou (1973) and Monahan (1984).
Each partial autocorrelation therefore takes an
[`rhobit_link`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html)
and every free vector gives a stationary, positive definite matrix.

The autocorrelations follow from the same recursion. Writing
\\\phi^{(k)}\\ for the coefficients of the order-\\k\\ predictor,
\$\$\phi^{(k)}\_k = r_k, \qquad \phi^{(k)}\_j = \phi^{(k-1)}\_j - r_k
\phi^{(k-1)}\_{k-j},\$\$ and the Yule-Walker equations at order \\k\\
give \\\rho_k = \sum\_{j\<k} \phi^{(k)}\_j \rho\_{k-j} + r_k\\, after
which \\\rho_h = \sum_j \phi_j \rho\_{h-j}\\ for every lag beyond the
order. The whole map from the partial autocorrelations to the matrix is
therefore **polynomial**, built from sums and products alone, and its
derivatives to fourth order are obtained by carrying jets through the
recursion rather than by expanding it.

Two quantities are closed form. The innovation variances of the
Levinson-Durbin recursion give \$\$\log\lvert M \rvert = p\log\gamma_0 +
\sum\_{k=1}^{q} (p - k)\log(1 - r_k^{2}),\$\$ one term per free value,
so the log-determinant is separable and every mixed derivative of it is
exactly zero. And the inverse is **banded of bandwidth \\q\\**: an
autoregression of order \\q\\ is Markov of that order, so its precision
carries no entry beyond the \\q\\-th off-diagonal. It is assembled from
the prediction form \\M^{-1} = U^\top D^{-1} U\\, with \\U\\ unit lower
triangular holding the predictor coefficients and \\D\\ the innovation
variances, rather than by a factorization.

[`ar1`](https://statmodels7.github.io/parameters7/reference/ar1.md) is
the case \\q = 1\\ written out: there the autocorrelation is simply
\\\rho^{h}\\, the determinant is \\(1-\rho^2)^{p-1}\\ and the inverse is
tridiagonal in three lines, so it keeps its own closed forms and does
not go through the recursion.

The name is `autoregressive()` rather than
[`ar()`](https://rdrr.io/r/stats/ar.html) because
[`ar`](https://rdrr.io/r/stats/ar.html) is a function of stats, and a
package meant to be attached alongside others should not mask one.

## References

Barndorff-Nielsen, O. and Schou, G. (1973). On the parametrization of
autoregressive models by partial autocorrelations. *Journal of
Multivariate Analysis* 3, 408-419.

Monahan, J. F. (1984). A note on enforcing stationarity in
autoregressive moving average models. *Biometrika* 71, 403-404.

## See also

[`ar1`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`compound_symmetry`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)

## Examples

``` r
s <- autoregressive(6, order = 2)
s@free_names
#> [1] "log_scale" "z_pacf1"   "z_pacf2"  

eta <- c(log(2), atanh(0.7), atanh(-0.3))
round(param_value(s, eta), 4)
#>         v1      v2     v3     v4      v5      v6
#> v1  2.0000  1.4000 0.6740 0.1933 -0.0263 -0.0819
#> v2  1.4000  2.0000 1.4000 0.6740  0.1933 -0.0263
#> v3  0.6740  1.4000 2.0000 1.4000  0.6740  0.1933
#> v4  0.1933  0.6740 1.4000 2.0000  1.4000  0.6740
#> v5 -0.0263  0.1933 0.6740 1.4000  2.0000  1.4000
#> v6 -0.0819 -0.0263 0.1933 0.6740  1.4000  2.0000

# the precision is banded of bandwidth two, the process being Markov of
# order two
round(param_solve(s, eta), 4)
#>         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
#> [1,]  1.0774 -0.9804  0.3232  0.0000  0.0000  0.0000
#> [2,] -0.9804  1.9695 -1.2745  0.3232  0.0000  0.0000
#> [3,]  0.3232 -1.2745  2.0665 -1.2745  0.3232  0.0000
#> [4,]  0.0000  0.3232 -1.2745  2.0665 -1.2745  0.3232
#> [5,]  0.0000  0.0000  0.3232 -1.2745  1.9695 -0.9804
#> [6,]  0.0000  0.0000  0.0000  0.3232 -0.9804  1.0774

# the round trip closes exactly
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 1.110223e-16
```
