# All Orderings of a Tuple

Returns the distinct permutations of an index tuple as a list. When the
tuple has repeated entries, permutations that coincide appear once, so
`combinat_perms(c(1, 1))` gives one element.
[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md)
restores the multiplicity with the product of the factorials of the
repeat counts.

## Usage

``` r
combinat_perms(x)
```

## Arguments

- x:

  An integer vector, of length 1 to 4 in every call the package makes.

## Value

A list of the distinct orderings of `x`, each an integer vector:
`factorial(length(x))` of them when the entries of `x` are all
different, and fewer when some repeat.

## See also

[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
the only caller.
