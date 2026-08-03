# The Number of Diagonal Entries Each Free Value Owns

One for an unshared structure, and the whole diagonal for a shared one.

## Usage

``` r
diag_multiplicity(s)
```

## Arguments

- s:

  A
  [`DiagStruct`](https://statmodels7.github.io/covstructs7/reference/DiagStruct.md)
  object.

## Value

An integer vector of length `s@n_free`.
