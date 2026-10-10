# Construct an Autoregressive Parameter

Returns an object holding the covariance of \\p\\ consecutive
observations of a stationary autoregression of order \\q\\,

\$\$y_t = \phi_1 y\_{t-1} + \cdots + \phi_q y\_{t-q} +
\varepsilon_t,\$\$

parametrized by its marginal variance \\\gamma_0\\ and its \\q\\ partial
autocorrelations \\r_1, \dots, r_q\\. That is \\q + 1\\ free values
**whatever the dimension**: `n_free` is 3 for an order-2 process,
whether it is observed 6 times or 200 times.

## Usage

``` r
autoregressive(dimension, order, link_scale = linkfunctions7::log_link())
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

## Value

An object of class
[`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md),
with `n_free` equal to \\q + 1\\, `free_names` the tagged `log_scale`,
`z_pacf1`, ..., `z_pacfq`, `rank` equal to `dimension`, an empty
`null_basis`, and `param_name` `"ar(<q>)"` with the order written in,
for example `"ar(2)"`.

## Why the partial autocorrelations carry the parametrization

The coefficients \\\phi\\ are stationary exactly when the roots of \\1 -
\phi_1 z - \cdots - \phi_q z^{q}\\ lie outside the unit circle, and that
set is not a box. At \\q = 2\\ it is the open triangle with vertices
\\(-2, -1)\\, \\(2, -1)\\ and \\(0, 1)\\, of area 4 inside a bounding
box \\(-2, 2) \times (-1, 1)\\ of area 8: **exactly half the box is
non-stationary**, and \\\phi = (1.5, 0.6)\\, which sits comfortably
inside the box, has a root of modulus 0.547. Scalar links onto intervals
therefore cannot cover the stationary region, for any choice of
intervals.

The partial autocorrelations do not have this problem. Each lies in
\\(-1, 1)\\ independently of the others, and the Levinson-Durbin
recursion carries them onto the stationary coefficients bijectively,
which is the transformation of Barndorff-Nielsen and Schou (1973) and
Monahan (1984). Each takes a
[`linkfunctions7::rhobit_link()`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html),
and every free vector then gives a stationary positive definite matrix.

In double precision that last statement holds only away from the edge of
the chart. As the partial autocorrelations approach \\\pm 1\\, the
smallest eigenvalue of the matrix falls to the rounding level, and with
several free values of absolute size 6 to 10, depending on \\p\\ and
\\q\\, the computed matrix can be indefinite. The log-determinant, its
derivatives and the innovation variances of the inverse use the factors
\\1 - r_k^2\\, which are evaluated from the free values by
[`sech2()`](https://statmodels7.github.io/parameters7/reference/sech2.md)
and
[`log_sech2()`](https://statmodels7.github.io/parameters7/reference/sech2.md)
and keep their accuracy in that region.

## The polynomial map from the partial autocorrelations

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
recursion in compiled code, the product rule written out per order, with
nothing differenced.

## The log-determinant and the inverse

The innovation variances of the Levinson-Durbin recursion give

\$\$\log\lvert M \rvert = p\log\gamma_0 + \sum\_{k=1}^{q} (p -
k)\log(1 - r_k^{2}),\$\$

one term per free value, so the log-determinant is **separable** and
every mixed derivative of it is exactly zero: at order 4 and \\q = 2\\,
12 of the 15 components are mixed and are 0 by construction. It costs a
sum of \\q + 1\\ terms at any dimension.

The inverse is **banded of bandwidth \\q\\**: an autoregression of order
\\q\\ is Markov of that order, so its precision carries no entry beyond
the \\q\\-th off-diagonal, and every entry outside the band is exactly
0. It is assembled from the prediction form \\M^{-1} = U^\top D^{-1}
U\\, with \\U\\ unit lower triangular holding the predictor coefficients
and \\D\\ the innovation variances, so no factorization is taken.

## Implementation

Each derivative order has its own compiled kernel, reached through
[`ar_tables()`](https://statmodels7.github.io/parameters7/reference/ar_tables.md),
which returns the components of that order only, and the Toeplitz
matrices are filled in compiled code.

## Against ar1()

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) is
the case \\q = 1\\ written out: there the autocorrelation is
\\\rho^{h}\\, the determinant is \\(1-\rho^2)^{p-1}\\ and the inverse is
tridiagonal in three lines, so it keeps its own closed forms and does
not go through the recursion. The two agree to rounding in the value,
the inverse, the log-determinant and the four derivative orders. The
free name of
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) is
`z_rho` where this one's is `z_pacf1`, and at \\q = 1\\ the two are the
same number.

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
which reports the autoregressive coefficients that this parametrization
does not carry.

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
# with the Jacobian that a delta method needs.
param_readable(s, eta)$value
#> scale pacf1 pacf2  phi1  phi2 
#>  2.00  0.70 -0.30  0.91 -0.30 
```
