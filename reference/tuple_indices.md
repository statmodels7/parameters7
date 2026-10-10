# The Index Tuples of a Given Width

The enumeration behind
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md),
taken over a plain count of variables instead of over a parameter
object, so that anything holding derivatives over \\d\\ variables can
use it without constructing a
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).
It forwards to
[`numericals7::tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.html),
so the two enumerations cannot disagree.

## Usage

``` r
tuple_indices(d, order = 2L)
```

## Arguments

- d:

  The number of variables, a single non-negative integer.

- order:

  The derivative order: 1, 2, 3 or 4. Any other value signals the error
  `'order' must be 1, 2, 3 or 4.`

## Value

A list of `choose(d + order - 1, order)` integer vectors, each of length
`order`, holding positions in `1:d` in non-decreasing order within a
tuple. At order 2 the diagonal pairs come first.

## See also

[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md),
the parameter-facing wrapper, and
[`numericals7::tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.html),
the implementation.
