# The Factor of a Log-Cholesky Parameter, and Its Derivatives

Returns \\\partial^S L\\ for a multiset \\S\\ of free-value indices, or
`NULL` when that derivative is identically zero.

## Usage

``` r
chol_dfactor(s, l, ks)
```

## Arguments

- s:

  A
  [`LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- l:

  The factor at the point, from
  [`chol_assemble`](https://statmodels7.github.io/parameters7/reference/chol_assemble.md).

- ks:

  A multiset of free-value indices, possibly empty; the empty one gives
  \\L\\ itself.

## Value

A numeric matrix, or `NULL`.

## Details

A free value below the diagonal enters \\L\\ linearly, so its second
derivative vanishes; a diagonal one enters through its exponential, so
every repeated derivative in that same value regenerates \\L\_{ii}
E\_{ii}\\. Everything else is zero, which is what makes the Leibniz sum
of
[`chol_leibniz`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md)
short.

## See also

[`chol_leibniz`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
[`leibniz_gram`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
