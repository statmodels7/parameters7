# The Eigendecomposition and Divided-Difference Tables at a Point

Everything the derivative contractions need: the eigendecomposition of
\\S\\, the rotated basis directions, and the divided-difference tables
of the orders asked for, computed once per free vector.

## Usage

``` r
mlog_tables(s, eta, order)
```

## Arguments

- s:

  A
  [`MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The highest derivative order wanted.

## Value

A list with `q`, `lam`, `e` (rotated directions) and the tables `dd2` to
`dd5` up to `order + 1` points.
