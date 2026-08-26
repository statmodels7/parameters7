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
  or `NULL` where that derivative is identically zero. The empty vector
  must give \\L\\ itself.

- tuple:

  The index tuple, an integer vector whose length is the derivative
  order. Repeated entries are handled correctly, the sum running over
  subsets of positions.

- p:

  The side of the matrix.

## Value

A symmetric `p` by `p` numeric matrix, with no dimnames.

## Details

The Leibniz rule distributes the differentiations of a product over its
two factors in every way, so \$\$\partial^T (L L^\top) = \sum\_{S
\subseteq T} (\partial^S L)(\partial^{T \setminus S} L)^\top,\$\$ the
sum running over subsets of *positions* in the tuple, which handles a
repeated index correctly without a multiplicity bookkeeping of its own.

The result is symmetric by construction, with no symmetrizing step: the
term for a subset \\S\\ and the term for its complement are transposes
of each other, so the sum pairs off. Measured at orders 1 to 3 on a
random factor, the asymmetry is exactly 0.

The loop walks the \\2^{\|T\|}\\ subsets through a bit mask, so the cost
is \\2^{\text{order}}\\ matrix products at worst, and far fewer in
practice: a term whose `dfactor` returns `NULL` is skipped before the
product is formed, and for both families that use this most terms do.

## See also

[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
for the other piece of shared arithmetic,
[`chol_dfactor()`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md)
and
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
for the two `dfactor` arguments the package supplies, and
[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
the compiled route that replaces this for the log-Cholesky family.
