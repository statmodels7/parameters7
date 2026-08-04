# Names of the Distinct Derivative Components

The keys of the derivative lists: one per unordered tuple of free
values. At order 2, the keys of
[`param_d2`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
and
[`param_d2logdet`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
diagonal pairs first; at orders 3 and 4, the keys of
[`param_d3`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
lexicographic combinations with repetition.

## Usage

``` r
param_tuple_names(s, order = 2L)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- order:

  The derivative order: 2 (default), 3 or 4.

## Value

A character vector, one entry per distinct component.

## Details

This exists so that nothing has to recover a tuple by splitting a key
apart. Generating the keys and the index tuples from one enumeration
cannot be fooled by a free name that contains the separator, which
splitting can. Use
[`param_tuple_indices`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
for the tuples themselves.

## See also

[`param_tuple_indices`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)

## Examples

``` r
param_tuple_names(log_cholesky(2))
#> [1] "log_L1:log_L1" "log_L2:log_L2" "L2.1:L2.1"     "log_L1:log_L2"
#> [5] "log_L1:L2.1"   "log_L2:L2.1"  
param_tuple_names(scalar_matrix(2), 3)
#> [1] "log_scale:log_scale:log_scale"
```
