# The Index Tuples of a Given Width

The enumeration behind
[`param_tuple_indices`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md),
taken over a number of variables rather than over a parameter, so that
anything holding derivatives over \\d\\ variables – a jet, for instance
– can share it.

## Usage

``` r
tuple_indices(d, order = 2L)
```

## Arguments

- d:

  The number of variables.

- order:

  The derivative order, 1 to 4.

## Value

A list of integer vectors of length `order`.

## See also

[`param_tuple_indices`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
