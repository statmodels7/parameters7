# Ordered Set Partitions of the First k Positions

Every way of splitting `1:k` into non-empty blocks **and ordering the
blocks**, which is what the derivative of an inverse sums over: the
factors are matrices and do not commute, so two orderings of the same
split are two terms.

## Usage

``` r
ordered_set_partitions(k)
```

## Arguments

- k:

  A single positive integer, at most 4 in this package's use.

## Value

A list of lists of integer vectors, each inner list one ordered
partition.

## Details

Built from
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html),
the toolkit's single copy of the unordered enumeration, by taking every
permutation of each partition's blocks. The count is the Fubini number:
1, 3, 13, 75 for `k` of 1 to 4.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md),
the only caller.
