# Names of the Distinct Second-Derivative Components

The keys of
[`struct_d2matrix`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.md)
and
[`struct_d2logdet`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.md):
one per unordered pair of free values, diagonal first and then the
off-diagonal pairs in lexicographic order.

## Usage

``` r
struct_pair_names(s)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

## Value

A character vector of length \\d(d+1)/2\\.

## Details

This exists so that nothing has to recover a pair by splitting a key
apart. Generating the keys and the index pairs from one enumeration
cannot be fooled by a free name that contains the separator, which
splitting can. Use
[`struct_pair_indices`](https://statmodels7.github.io/covstructs7/reference/struct_pair_indices.md)
for the pairs themselves.

## See also

[`struct_pair_indices`](https://statmodels7.github.io/covstructs7/reference/struct_pair_indices.md)

## Examples

``` r
struct_pair_names(log_cholesky(2))
#> [1] "log_L1:log_L1" "log_L2:log_L2" "L2.1:L2.1"     "log_L1:log_L2"
#> [5] "log_L1:L2.1"   "log_L2:L2.1"  
```
