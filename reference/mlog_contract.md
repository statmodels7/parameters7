# One Derivative Component by the Daleckii-Krein Contraction

Contracts the rotated directions of one index tuple against the
divided-difference table of the matching order, summing over the
orderings of the directions, and rotates back.

## Usage

``` r
mlog_contract(tb, dirs)
```

## Arguments

- tb:

  The tables of
  [`mlog_tables`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md).

- dirs:

  The free-value indices of the tuple, possibly repeated.

## Value

A symmetric numeric matrix.
