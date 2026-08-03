# The Diagonal Entries Behind a Free Vector

Applies the structure's link, recycling a shared value across the
diagonal.

## Usage

``` r
diag_entries(s, eta)
```

## Arguments

- s:

  A
  [`DiagStruct`](https://statmodels7.github.io/covstructs7/reference/DiagStruct.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A numeric vector of length `s@dimension`.
