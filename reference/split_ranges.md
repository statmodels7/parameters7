# Consecutive Index Ranges of Given Widths

Turns the widths \\n_1, \ldots, n_B\\ into the ranges they occupy when
laid end to end, as a list of integer vectors. A width of zero gives an
empty range rather than being dropped, so that the list stays aligned
with the blocks.

## Usage

``` r
split_ranges(widths)
```

## Arguments

- widths:

  An integer vector of widths.

## Value

A list of integer vectors, one per width.
