# Derivative Components of a Log-Cholesky Parameter

Assembles one derivative order by the Leibniz rule on \\M = LL^ op\\,
the factor's derivatives coming from
[`chol_dfactor`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md).

## Usage

``` r
chol_leibniz(s, eta, order)
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

## See also

[`leibniz_gram`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
