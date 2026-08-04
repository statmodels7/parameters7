# The Factor of a Correlation Parameter, and Its Derivatives

Returns \\\partial^S L\\ for a multiset \\S\\ of free-value indices, or
`NULL` when that derivative is identically zero.

## Usage

``` r
corr_dfactor(s, tb, ks)
```

## Arguments

- s:

  A
  [`CorrelationParam`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- tb:

  The tables of
  [`corr_tables`](https://statmodels7.github.io/parameters7/reference/corr_tables.md).

- ks:

  A multiset of free-value indices, possibly empty.

## Value

A numeric matrix, or `NULL`.

## Details

An entry of \\L\\ depends only on the angles of its own row, so a
multiset spanning two rows gives zero; within a row, an entry gives zero
unless every differentiated angle appears among its factors.
