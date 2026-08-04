# Index Tuples Behind the Derivative Component Names

The unordered index tuples
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
names, in exactly the same order. At order 2 the diagonal pairs come
first and then the off-diagonal ones, because consumers index a Hessian
that way; at orders 3 and 4 the tuples are the lexicographic
combinations with repetition, matching the enumeration distributions7
uses for its higher derivatives.

## Usage

``` r
param_tuple_indices(s, order = 2L)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- order:

  The derivative order: 1 to 4.

## Value

A list of integer vectors of length `order`.

## See also

[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)

## Examples

``` r
param_tuple_indices(log_cholesky(2), 2)
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
length(param_tuple_indices(log_cholesky(2), 3))
#> [1] 10
```
