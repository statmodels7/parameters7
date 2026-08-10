# Orderings of a Multiset, Counted With Multiplicity

All \\n!\\ orderings of a vector of indices, without deduplicating the
ones that coincide because an index repeats.

## Usage

``` r
multiset_orderings(v)
```

## Arguments

- v:

  An integer vector.

## Value

A list of integer vectors.

## Details

The distinction is load bearing. The cyclic sum behind the
log-determinant expansion runs over \\(n-1)!\\ orderings, and two that
happen to be equal still count twice; deduplicating them makes the third
derivative of a component in one weight too small by a factor of two and
the fourth by six.
