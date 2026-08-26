# Consecutive Index Ranges of Given Widths

Turns the widths \\n_1, \ldots, n_B\\ into the ranges they occupy when
laid end to end. It is called twice at construction, once for the
blocks' rows in the matrix and once for their stretches of the free
vector.

## Usage

``` r
split_ranges(widths)
```

## Arguments

- widths:

  An integer vector of non-negative widths.

## Value

A list of integer vectors, one per width, in order.

## Details

A width of zero gives an **empty range**, `integer(0)`, instead of being
dropped, so the list stays aligned with the blocks and the \\b\\-th
element is always the \\b\\-th block's. That case is reachable: a fully
known
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
has no free values at all, so `split_ranges(c(2, 0, 3))` returns `1:2`,
`integer(0)` and `3:5`.

## See also

[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
the only caller.
