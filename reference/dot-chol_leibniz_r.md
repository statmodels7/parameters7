# The R Twin of the Compiled Leibniz Assembly

The same components as `chol_leibniz_cpp`, built through the dense
products of
[`leibniz_gram`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md).
It is not the production route: the compiled kernel exploits that every
derivative of the factor is a single-entry matrix, so each Leibniz term
is one row, one column or one cell, and at \\p = 8\\, order 4, it
measured four orders of magnitude faster. The twin is kept as the
independent reference the tests compare against.

## Usage

``` r
.chol_leibniz_r(s, eta, order)
```

## Arguments

- s:

  A
  [`LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order.

## Value

A named list of symmetric matrices.
