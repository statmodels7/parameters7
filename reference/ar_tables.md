# The Levinson-Durbin Recursion and One Derivative Order

Runs the compiled kernel of one order: `ar_value_cpp()` at order 0 and
`ar_d1_cpp()` to `ar_d4_cpp()` above. The scale and the partial
autocorrelations enter as their link inverses with derivatives to
`order`, and the autocovariances and the autoregressive coefficients
come back with the derivative components of that order alone.

## Usage

``` r
ar_tables(s, eta, order)
```

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order, an integer from 0 to 4.

## Value

A list with `gamma` and `phi`. At order 0, `gamma` holds the
autocovariances at lags \\0, \dots, p-1\\ and `phi` the \\q\\
coefficients. Above, each is a matrix (\\p\\ and \\q\\ rows) with one
column per index tuple of
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
at that order, in its order.

## Details

The coefficients are multilinear in the partial autocorrelations \\r\\:
step \\k\\ of the recursion is affine in \\r_k\\ and the earlier
coefficients do not involve it. The kernel differentiates the recursion
in \\r\\ for every multiset of indices up to `order`, the products by
Leibniz over the sub-multisets without a repeated index, and applies the
links once at the end through the partial Bell polynomials, each \\r_k\\
depending on its own free value alone.

## See also

[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the recursion itself, and
[`ar_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_derivative.md),
which assembles the matrices.
