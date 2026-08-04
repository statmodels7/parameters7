# The Diagonal Entries Behind a Free Vector

Applies the parameter's link, recycling a shared value across the
diagonal.

## Usage

``` r
diag_entries(s, eta)
```

## Arguments

- s:

  A
  [`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A numeric vector of length `s@dimension`.
