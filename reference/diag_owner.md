# The Free Value Each Diagonal Entry Belongs To

The index into the free vector of the value controlling each entry: the
identity for an unshared parameter, and all ones for a shared one.

## Usage

``` r
diag_owner(s)
```

## Arguments

- s:

  A
  [`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

## Value

An integer vector of length `s@dimension`.
