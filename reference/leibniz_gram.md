# A Gram Product's Derivatives From Its Factor's

The derivative of \\M = L L^\top\\ for one index tuple, from the
derivatives of \\L\\.

## Usage

``` r
leibniz_gram(dfactor, tuple, p)
```

## Arguments

- dfactor:

  A function of a (possibly empty, possibly repeating) integer vector of
  free-value indices, returning the corresponding derivative of \\L\\,
  or `NULL` when that derivative is identically zero. The empty vector
  must give \\L\\ itself.

- tuple:

  The index tuple.

- p:

  The side of the matrix.

## Value

A symmetric numeric matrix.

## Details

The Leibniz rule distributes the differentiations of a product over its
two factors in every way, so \$\$\partial^T (L L^\top) = \sum\_{S
\subseteq T} (\partial^S L)(\partial^{T \setminus S} L)^\top,\$\$ the
sum running over subsets of *positions* in the tuple, which handles a
repeated index correctly without a multiplicity bookkeeping of its own.

## See also

[`compose4`](https://statmodels7.github.io/parameters7/reference/compose4.md)
