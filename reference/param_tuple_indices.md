# Index Tuples Behind the Derivative Component Names

Returns the unordered index tuples that
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
labels, in exactly the same order, as integer positions into the free
vector. A consumer that has to know *which* free values a component
differentiates in reads these; a consumer that only has to display a
label reads the names. At order 2 the diagonal pairs come first and the
off-diagonal ones after, which is how a Hessian is filled; at orders 3
and 4 the tuples are the lexicographic combinations with repetition.

## Usage

``` r
param_tuple_indices(s, order = 2L)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md),
  whose `n_free` supplies \\d\\.

- order:

  The derivative order: `1`, `2` (the default), `3` or `4`. Anything
  else throws `'order' must be 1, 2, 3 or 4.`

## Value

A list of `choose(s@n_free + order - 1, order)` integer vectors, each of
length `order`, holding positions in `1:s@n_free` in non-decreasing
order within a tuple.

## Details

The order at orders 3 and 4 matches the enumeration distributions7 uses
for its own higher derivatives, so a consumer contracting a parameter's
derivative array against a distribution's can walk the two lists
together without reindexing. Both come from
[`numericals7::tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.html),
which is the single copy of the enumeration in the toolkit.

## Notation

\\d\\ is the length of the free vector, `s@n_free`, and \\k\\ the
derivative order.

## See also

[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
for the same tuples as labels, and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the lists they index.

## Examples

``` r
s <- log_cholesky(2)

# Order 2: the three diagonal pairs first, then the three off-diagonal ones.
param_tuple_indices(s, 2)
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2 2
#> 
#> [[3]]
#> [1] 3 3
#> 
#> [[4]]
#> [1] 1 2
#> 
#> [[5]]
#> [1] 1 3
#> 
#> [[6]]
#> [1] 2 3
#> 

# The tuples and the names are the same enumeration, so they line up.
idx <- param_tuple_indices(s, 3)
nm <- param_tuple_names(s, 3)
identical(nm, vapply(idx, function(t) paste(s@free_names[t], collapse = ":"),
                     character(1)))
#> [1] TRUE

# What they are for: reading a component's indices to contract with.
d2 <- param_d2(s, c(0.2, -0.1, 0.4))
d1 <- param_d1(s, c(0.2, -0.1, 0.4))
i <- which(param_tuple_names(s, 2) == "log_L1:L2.1")
param_tuple_indices(s, 2)[[i]]
#> [1] 1 3
dim(d2[[i]])
#> [1] 2 2

# Counts, over unordered tuples.
vapply(1:4, function(k) length(param_tuple_indices(s, k)), integer(1))
#> [1]  3  6 10 15
```
