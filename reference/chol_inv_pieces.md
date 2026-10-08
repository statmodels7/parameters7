# The Pieces Behind a Log-Cholesky Inverse's Derivatives

\\L\\, \\G = L^{-1}\\, and for each free value \\B_k = G\\\partial_k L\\
and \\C_k = B_k + B_k^\top\\, which
[`param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
and
[`param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
read for a
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
parameter.

## Usage

``` r
chol_inv_pieces(s, eta)
```

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  Its free vector.

## Value

A list with `l`, `g`, `b` and `c`, the last two lists over the free
values.
