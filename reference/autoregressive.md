# Construct an Autoregressive Parameter

Returns an object holding the covariance of \\p\\ consecutive
observations of a stationary autoregression of order \\q\\,

\$\$y_t = \phi_1 y\_{t-1} + \cdots + \phi_q y\_{t-q} +
\varepsilon_t,\$\$

parametrized by its marginal variance \\\gamma_0\\ and its \\q\\ partial
autocorrelations \\r_1, \dots, r_q\\. That is \\q + 1\\ free values
**whatever the dimension**: measured, `n_free` is 3 for an order-2
process observed 6 times and 3 for one observed 200 times.

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
  Must exceed `order`, since a stretch of \\q\\ observations does not
  identify \\q\\ partial autocorrelations; \\p = q + 1\\ is legal and is
  the smallest legal dimension.

- order:

  The order \\q\\ of the autoregression: a single positive whole number,
  no smaller than 1.

- link_scale:

  A linkfunctions7 link carrying the marginal variance,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map onto the positive half line, so
  `identity_link()` is rejected, and from the whole real line, which
  rules out `sqrt_link()` and its relatives; see
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  for the two conditions.

- role:

  A label recording which side of a model the matrix parametrizes:
  `"either"` (the default), `"covariance"` or `"precision"`. No numeric
  result depends on it.

## Value

An object of class
[`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md),
with `n_free` equal to \\q + 1\\, `free_names` the tagged `log_scale`,
`z_pacf1`, ..., `z_pacfq`, `rank` equal to `dimension`, an empty
`null_basis`, and `param_name` `"ar(q)"`.

## Why the partial autocorrelations carry the parametrization

The coefficients \\\phi\\ are stationary exactly when the roots of \\1 -
\phi_1 z - \cdots - \phi_q z^{q}\\ lie outside the unit circle, and that
set is not a box. At \\q = 2\\ it is the open triangle with vertices
\\(-2, -1)\\, \\(2, -1)\\ and \\(0, 1)\\, of area 4 inside a bounding
box \\(-2, 2) \times (-1, 1)\\ of area 8: **exactly half the box is
non-stationary**, and \\\phi = (1.5, 0.6)\\, which sits comfortably
inside the box, has a root of modulus 0.547. So no collection of scalar
links onto intervals can cover the region, whatever intervals are
chosen.

The partial autocorrelations do not have this problem. Each lies in
\\(-1, 1)\\ independently of the others, and the Levinson-Durbin
recursion carries them onto the stationary coefficients bijectively,
which is the transformation of Barndorff-Nielsen and Schou (1973) and
Monahan (1984). Each takes a
[`linkfunctions7::rhobit_link()`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html),
and every free vector then gives a stationary positive definite matrix.

In double precision that last statement has a boundary, and it belongs
to the chart, not to this family. At \\\lvert \eta_k \rvert = 6\\ the
matrix is strictly positive definite with an eigenvalue ratio of \\5
\times 10^{-12}\\; by \\\lvert \eta_k \rvert = 10\\ the partial
autocorrelation is \\1 - 4 \times 10^{-9}\\ and the smallest eigenvalue
has crossed zero at the rounding floor, \\-3 \times 10^{-17}\\ of the
largest.

## The map is polynomial, so nothing is differenced

The autocorrelations follow from the same recursion. Writing
\\\phi^{(k)}\\ for the coefficients of the order-\\k\\ predictor,

\$\$\phi^{(k)}\_k = r_k, \qquad \phi^{(k)}\_j = \phi^{(k-1)}\_j - r_k
\phi^{(k-1)}\_{k-j},\$\$

and the Yule-Walker equations at order \\k\\ give \\\rho_k =
\sum\_{j\<k} \phi^{(k)}\_j \rho\_{k-j} + r_k\\, after which \\\rho_h =
\sum_j \phi_j \rho\_{h-j}\\ for every lag beyond the order. The whole
map from the partial autocorrelations to the matrix is therefore
**polynomial**, built from sums and products alone, and its derivatives
to fourth order come from propagating the derivative arrays through the
recursion in compiled code, the product rule written out per order.
Measured against one central difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
the first derivatives agree to \\5 \times 10^{-11}\\, which is the
difference's own accuracy.

## Two quantities are closed form

The innovation variances of the Levinson-Durbin recursion give

\$\$\log\lvert M \rvert = p\log\gamma_0 + \sum\_{k=1}^{q} (p -
k)\log(1 - r_k^{2}),\$\$

one term per free value, so the log-determinant is **separable** and
every mixed derivative of it is exactly zero: at order 4 and \\q = 2\\,
12 of the 15 components are 0 by construction. It agrees with
[`determinant()`](https://rdrr.io/r/base/det.html) to \\2 \times
10^{-15}\\ at \\p = 5\\ and \\4 \times 10^{-14}\\ at \\p = 100\\, and
costs a sum of \\q + 1\\ terms at either size.

The inverse is **banded of bandwidth \\q\\**: an autoregression of order
\\q\\ is Markov of that order, so its precision carries no entry beyond
the \\q\\-th off-diagonal. Measured at \\q = 3\\ and \\p = 9\\, every
one of the 30 entries outside the band is exactly 0, not merely small.
It is assembled from the prediction form \\M^{-1} = U^\top D^{-1} U\\,
with \\U\\ unit lower triangular holding the predictor coefficients and
\\D\\ the innovation variances, so no factorization is taken.

## What a call costs

[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
runs the recursion once and packs **every** order up to the fourth, so
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
pays most of what
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
pays and the difference between them is the R-level assembly. Seconds
per call, over repetition loops sized by elapsed time:

|       |       |               |            |            |
|-------|-------|---------------|------------|------------|
| \\q\\ | \\p\\ | `param_value` | `param_d1` | `param_d4` |
| 1     | 10    | 0.00018       | 0.00024    | 0.00033    |
| 1     | 200   | 0.00143       | 0.00191    | 0.00375    |
| 2     | 200   | 0.00398       | 0.00516    | 0.01203    |
| 4     | 200   | 0.03313       | 0.03563    | 0.08313    |

## Against ar1()

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) is
the case \\q = 1\\ written out: there the autocorrelation is
\\\rho^{h}\\, the determinant is \\(1-\rho^2)^{p-1}\\ and the inverse is
tridiagonal in three lines, so it keeps its own closed forms and does
not go through the recursion. The two agree: at \\p = 8\\ and \\\rho =
0.6\\ the values differ by \\1.4 \times 10^{-17}\\, the log-determinants
by 0, the inverses by \\2.2 \times 10^{-16}\\ and the four derivative
orders by \\1.1 \times 10^{-16}\\ to \\1.1 \times 10^{-14}\\. Cost is
not the reason to prefer one: they measure 0.00033 s and 0.00050 s for
one fourth-order array.
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)'s
free name is `z_rho` where this one's is `z_pacf1`, and at \\q = 1\\ the
two are the same number.

## The name

The name is `autoregressive()` rather than
[`ar()`](https://rdrr.io/r/stats/ar.html) because
[`stats::ar()`](https://rdrr.io/r/stats/ar.html) is a function of stats,
and a package meant to be attached alongside others should not mask one.

## Notation

\\p\\ is the dimension, \\q\\ the order, \\\gamma_0\\ the marginal
variance, \\r_k\\ the \\k\\-th partial autocorrelation, \\\phi_j\\ the
autoregressive coefficients, \\\rho_h\\ the autocorrelation at lag
\\h\\, and \\\eta = (\log \gamma_0, \mathrm{atanh}\\ r_1, \dots)\\ the
free vector.

## References

Barndorff-Nielsen, O. and Schou, G. (1973). On the parametrization of
autoregressive models by partial autocorrelations. *Journal of
Multivariate Analysis* **3**, 408-419.

Monahan, J. F. (1984). A note on enforcing stationarity in
autoregressive moving average models. *Biometrika* **71**, 403-404.

## See also

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
for \\q = 1\\ written out,
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
for the other two-value family, and
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md),
which reports the autoregressive coefficients this parametrization does
not carry.

## Examples

``` r
s <- autoregressive(6, order = 2)
s@free_names
#> [1] "log_scale" "z_pacf1"   "z_pacf2"  

eta <- c(log(2), atanh(0.7), atanh(-0.3))
M <- param_value(s, eta)
round(M, 4)
#>         v1      v2     v3     v4      v5      v6
#> v1  2.0000  1.4000 0.6740 0.1933 -0.0263 -0.0819
#> v2  1.4000  2.0000 1.4000 0.6740  0.1933 -0.0263
#> v3  0.6740  1.4000 2.0000 1.4000  0.6740  0.1933
#> v4  0.1933  0.6740 1.4000 2.0000  1.4000  0.6740
#> v5 -0.0263  0.1933 0.6740 1.4000  2.0000  1.4000
#> v6 -0.0819 -0.0263 0.1933 0.6740  1.4000  2.0000

# The first partial autocorrelation is the lag-one correlation.
c(from_matrix = M[1, 2] / M[1, 1], from_eta = tanh(eta[2]))
#> from_matrix    from_eta 
#>         0.7         0.7 

# The precision is banded of bandwidth two, the process being Markov of
# order two: everything beyond the second off-diagonal is exactly zero.
lag <- abs(outer(1:6, 1:6, "-"))
max(abs(param_solve(s, eta)[lag > 2]))
#> [1] 0

# The log-determinant is a sum of three terms, whatever p is.
c(closed = param_logdet(s, eta),
  from_determinant = as.numeric(determinant(M, logarithm = TRUE)$modulus))
#>           closed from_determinant 
#>        0.4149176        0.4149176 

# And it is separable, so every mixed derivative of it is exactly zero.
tuples <- param_tuple_indices(s, 2L)
mixed <- vapply(tuples, function(t) t[1] != t[2], logical(1))
max(abs(param_d2logdet(s, eta)[mixed]))
#> [1] 0

# The round trip closes.
max(abs(param_free(s, M) - eta))
#> [1] 1.110223e-16

# The coefficients are not free values, and param_readable() reports them
# with the Jacobian a delta method needs.
param_readable(s, eta)$value
#> scale pacf1 pacf2  phi1  phi2 
#>  2.00  0.70 -0.30  0.91 -0.30 
```
