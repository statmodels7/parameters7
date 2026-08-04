# The Number of Diagonal Entries Each Free Value Owns

One for an unshared parameter, and the whole diagonal for a shared one.

## Usage

``` r
diag_multiplicity(s)
```

## Arguments

- s:

  A
  [`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

## Value

An integer vector of length `s@n_free`.
