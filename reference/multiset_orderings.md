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

The cyclic sum behind the log-determinant expansion runs over \\(n-1)!\\
orderings, and two that happen to be equal still count twice. Counting
them once would make the trace term whose indices are all equal too
small by the factors \\2!\\ and \\3!\\, the numbers of orderings of the
tail at third and fourth order.

## See also

[`sum_struct_trace_term()`](https://statmodels7.github.io/parameters7/reference/sum_struct_trace_term.md),
the only caller, and
[`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md),
which does the same for the matrix logarithm's contraction.
