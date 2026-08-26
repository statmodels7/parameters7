# Orderings of a Multiset, Counted With Multiplicity

All \\n!\\ orderings of a vector of indices, **without** deduplicating
the ones that coincide because an index repeats:
`multiset_orderings(c(1, 1))` returns two elements.

## Usage

``` r
multiset_orderings(v)
```

## Arguments

- v:

  An integer vector, of length 0 to 3 in every call the package makes,
  the expansion at order \\n\\ permuting the \\n-1\\ indices after the
  first.

## Value

A list of `factorial(length(v))` integer vectors.

## Details

The distinction is load bearing, and its cost is measured. The cyclic
sum behind the log-determinant expansion runs over \\(n-1)!\\ orderings,
and two that happen to be equal still count twice; deduplicating them
leaves a third derivative in one weight too small by exactly 2 and a
fourth by exactly 6, which is \\2!\\ and \\3!\\, the orderings of the
tail. Both are numbers a reader would accept without noticing.

## See also

[`sum_struct_trace_term()`](https://statmodels7.github.io/parameters7/reference/sum_struct_trace_term.md),
the only caller, and
[`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md),
which does the same for the matrix logarithm's contraction.
