# Construct the Precision of an AR(1) Process

Returns the family whose value is \\\Sigma(\eta)^{-1}\\ for \\\Sigma\\
an [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
matrix: the **tridiagonal** precision of a first-order autoregression.
It carries
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)'s
free vector unchanged, so its two coordinates are the scale and the
correlation of the process whose precision this is.

## Usage

``` r
ar1_inv(dimension, link_scale = linkfunctions7::log_link())
```

## Arguments

- dimension:

  The side of the matrix, at least 2.

- link_scale:

  The link carrying the AR(1) process's variance, passed to
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md).
  Must map onto the positive half line.

## Value

An object of class
[`Ar1InvParam()`](https://statmodels7.github.io/parameters7/reference/Ar1InvParam.md),
with `n_free` 2, `free_names`
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)'s
and `rank` equal to `dimension`.

## Which side is which

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) is
the AR(1) **pattern**, \\\sigma^2\rho^{\|i-j\|}\\. A family that takes a
matrix parameter decides what it does with it, so
`mvgaussian1_distrib(p, ar1(p))` is the process whose covariance is
AR(1) and `mvgaussian2_distrib(p, ar1(p))` the model whose precision has
those entries, which is a different and legitimate model. This family is
the third object in that picture: the matrix whose inverse is an AR(1),
so `mvgaussian2_distrib(p, ar1_inv(p))` is again the AR(1) process,
written on the precision side.

## The value

A first-order autoregression is Markov, so its precision carries no
entry beyond the first off-diagonal. With \\\tau = \sigma^{-2}\\ and \\w
= (1-\rho^2)^{-1}\\,

\$\$\Omega = \tau\\G(\rho), \qquad G\_{11} = G\_{pp} = w, \quad G\_{ii}
= (1+\rho^2)w, \quad G\_{i,i\pm 1} = -\rho w,\$\$

every other entry exactly zero. The value is read from
[`param_solve.Ar1Param()`](https://statmodels7.github.io/parameters7/reference/param_solve.Ar1Param.md),
which writes that matrix down and performs no factorization, so there is
one copy of it in the package.

## The derivatives

The value is a **product** of a function of the first free value and a
matrix function of the second, so a component with \\a\\ scale indices
and \\b\\ correlation indices is \\(\partial^a\tau)(\partial^b G)\\ and
no mixed expansion is needed. The scale factor is the reciprocal of the
link's inverse, chained through \\x\mapsto x^{-1}\\; the pattern's three
distinct entries are written in \\w\\ alone,

\$\$w^{(k)} = \frac{k!}{2}\left\\(1-\rho)^{-(k+1)} + (-1)^k
(1+\rho)^{-(k+1)}\right\\,\$\$

with \\G\_{11}^{(k)} = w^{(k)}\\, \\G\_{ii} = 2w - 1\\ so
\\G\_{ii}^{(k)} = 2w^{(k)}\\ above order zero, and \\G\_{i,i\pm1}^{(k)}
= -(\rho\\w^{(k)} + k\\w^{(k-1)})\\ by the Leibniz rule. Each is then
chained onto the free value through the correlation's link.

This is what the family adds over `inverse_of(ar1(p))`, which reaches
the same numbers through the ordered-block-partition sum: there a
fourth-order component is 75 products of five matrices, here it is one
elementwise product. The two agree to machine precision, which is what a
test asserts.

## The log-determinant

\\\log\lvert\Omega\rvert = -\log\lvert\Sigma\rvert\\, so the value and
all four orders are
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)'s
negated, and the closed form \\p\log\sigma^2 + (p-1)\log(1-\rho^2)\\ is
what is negated.

## Notation

\\p\\ is the matrix side, \\\sigma^2\\ the AR(1) process's variance,
\\\rho\\ its lag-one correlation, \\\tau = \sigma^{-2}\\ and \\w =
(1-\rho^2)^{-1}\\.

## See also

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
for the family this inverts,
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
for the order-\\q\\ case, and
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for any other family's inverse.

## Examples

``` r
s <- ar1_inv(5)
eta <- c(log(2), atanh(0.6))

# The value is tridiagonal: an AR(1) is Markov, so the precision carries no
# entry beyond the first off-diagonal.
round(param_value(s, eta), 4)
#>         v1      v2      v3      v4      v5
#> v1  0.7812 -0.4688  0.0000  0.0000  0.0000
#> v2 -0.4688  1.0625 -0.4688  0.0000  0.0000
#> v3  0.0000 -0.4688  1.0625 -0.4688  0.0000
#> v4  0.0000  0.0000 -0.4688  1.0625 -0.4688
#> v5  0.0000  0.0000  0.0000 -0.4688  0.7812

# And it is exactly the inverse of the AR(1) it names.
max(abs(param_value(s, eta) %*% param_value(ar1(5), eta) - diag(5)))
#> [1] 4.440892e-16

# The written-out derivatives agree with the general composition's.
g <- inverse_of(ar1(5))
max(abs(unlist(param_d4(s, eta)) - unlist(param_d4(g, eta))))
#> [1] 1.918465e-13
```
