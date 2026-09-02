# Construct the Precision of an Autoregression of Order q

Returns the family whose value is \\\Sigma(\eta)^{-1}\\ for \\\Sigma\\
an
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
matrix: the precision of an autoregression of order \\q\\, **banded of
bandwidth \\q\\**. It carries
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)'s
free vector unchanged, so its coordinates are the scale and the partial
autocorrelations of the process whose precision this is.

## Usage

``` r
autoregressive_inv(dimension, order, ...)
```

## Arguments

- dimension:

  The side of the matrix.

- order:

  The autoregressive order \\q\\, passed to
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md).

- ...:

  Further arguments passed to
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
  such as `link_scale`.

## Value

An object of class
[`AutoregressiveInvParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveInvParam.md),
with `n_free`, `free_names` and `dimension`
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)'s.

## The value and the log-determinant are closed

An autoregression of order \\q\\ is Markov of that order, so its
precision carries no entry beyond the \\q\\-th diagonal. The value comes
from
[`param_solve.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_solve.AutoregressiveParam.md),
which reads it off the prediction form of the process rather than
factorizing, and the log-determinant is
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)'s
negated. Both are therefore \\O(p)\\ in the entries that matter and
exact.

## The derivative arrays

Written out, from the same prediction form the value comes from:

\$\$\Omega = U^\top \mathrm{diag}(\tau)\\ U,\$\$

with \\U\\ unit lower triangular of bandwidth \\q\\, row \\t\\ carrying
the coefficients of the best linear predictor of \\y_t\\ from its
predecessors, and \\\tau_t = 1/v_t\\ the reciprocal innovation
variances. Two facts make every order exact without a new recursion.

The lower-order rows of \\U\\ are the coefficients of the SAME family at
that order, measured to 0, so their derivative arrays come from the
compiled Levinson-Durbin recursion of
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
run once per order, and a component differentiating in a partial
autocorrelation an order does not reach is exactly zero. And \\\tau_t\\
is a PRODUCT of one factor per free value, \\1/v_0\\ from the scale and
\\(1-r_j^2)^{-1}\\ from each correlation the prediction has reached, so
a mixed derivative of it is a product of univariate derivatives and is
exactly zero where it differentiates in a factor a row does not carry.

What is left is the Leibniz rule over three factors, taken twice so that
a component costs \\2^m\\ matrix products rather than \\3^m\\.

Measured against `inverse_of(autoregressive(p, q))`, which reaches the
same numbers through the ordered-block-partition sum, at order four:
**11.2x** at \\p = 6, q = 2\\, 11.4x at \\p = 20, q = 3\\ and **22.8x**
at \\p = 100\\, with the two agreeing to 7e-11 over six shapes and two
free vectors each. That agreement is what licenses the written-out
route, the two sharing no arithmetic.

## See also

[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the family this inverts,
[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
for the order-one case, and
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the general composition both of them specialize.

## Examples

``` r
s <- autoregressive_inv(6, order = 2)
eta <- c(log(1.5), atanh(0.5), atanh(-0.2))

# Banded of bandwidth two: the process is Markov of order two.
lag <- abs(outer(1:6, 1:6, "-"))
max(abs(param_value(s, eta)[lag > 2]))
#> [1] 0

# And it is the inverse of the autoregression it names.
max(abs(param_value(s, eta) %*% param_value(autoregressive(6, 2), eta) -
        diag(6)))
#> [1] 2.356654e-16
```
