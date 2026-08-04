# Higher Derivatives of a Scaled Log-Pseudo-Determinant

The derivative components of orders two to four of
\\r\log\lambda(\eta) + \log\lvert P\rvert\_+\\ in the single free value.

## Usage

``` r
scaled_dlog(s, eta, order)
```

## Arguments

- s:

  A
  [`ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- eta:

  A numeric vector of one free value.

- order:

  The derivative order, 2 to 4.

## Value

A named numeric vector of length one, keyed by
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, order)`.

## Details

The constant contributes nothing beyond order zero, so every component
is the rank times the corresponding derivative of \\\log\lambda\\, which
is Faa di Bruno's chain for the logarithm as in
[`diag_dlog`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md).
With one free value there is exactly one component per order.

## See also

[`scaled_matrix`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
[`diag_dlog`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
