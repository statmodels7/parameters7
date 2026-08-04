# The Cholesky Factor Behind a Free Vector

Assembles \\L\\ from the free vector: the diagonal is the exponential of
the first `dimension` values, the rest are placed below it.

## Usage

``` r
chol_assemble(s, eta)
```

## Arguments

- s:

  A
  [`LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A lower triangular numeric matrix.
