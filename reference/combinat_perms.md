# All Orderings of a Tuple

Returns every permutation of an index tuple as a list, including the
ones that coincide when the tuple has repeats: `combinat_perms(c(1, 1))`
gives two elements and never one. The multiplicity is deliberate: it is
what leaves
[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md)'s
sum over orderings correct for a repeated tuple without a correction
factor.

## Usage

``` r
combinat_perms(x)
```

## Arguments

- x:

  An integer vector, of length 1 to 4 in every call the package makes.

## Value

A list of `factorial(length(x))` integer vectors.

## See also

[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
the only caller.
