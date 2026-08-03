# Index Pairs Behind the Second-Derivative Names

The unordered index pairs
[`struct_pair_names`](https://statmodels7.github.io/covstructs7/reference/struct_pair_names.md)
names, in exactly the same order: the diagonal \\(k, k)\\ first, then
the off-diagonal pairs \\(k, l)\\ with \\k \< l\\ in lexicographic
order.

## Usage

``` r
struct_pair_indices(s)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

## Value

A list of integer vectors of length 2.

## See also

[`struct_pair_names`](https://statmodels7.github.io/covstructs7/reference/struct_pair_names.md)

## Examples

``` r
struct_pair_indices(log_cholesky(2))
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
```
