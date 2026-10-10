# One Derivative Component by the Daleckii-Krein Contraction

Contracts the rotated directions of one index tuple against the
divided-difference table of the matching order and sums over the
orderings of the directions. The result is in the eigenbasis of \\S\\;
the caller rotates it back by \\Q\\.

## Usage

``` r
mlog_contract(tb, dirs)
```

## Arguments

- tb:

  The tables of
  [`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md),
  carrying at least `length(dirs) + 1` points.

- dirs:

  The free-value indices of the tuple, possibly repeated. Its length is
  the derivative order.

## Value

A symmetric numeric matrix with the side of `tb$q`, in the eigenbasis of
\\S\\ and with no dimnames.

## Details

The multilinear form sums over all orderings of the directions. The loop
runs over the distinct orderings, which
[`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md)
returns, and the sum is multiplied by \\\prod_j m_j!\\, where the
\\m_j\\ are the multiplicities of the repeated indices, because each
distinct ordering occurs that many times among all of them. Without that
factor the result would be too small by exactly that factor.

The cost is the reason this family's higher orders are expensive: at
order 4 a component with four different indices has 24 orderings, each a
contraction of \\O(p^5)\\ operations.

## See also

[`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md)
for the input,
[`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md)
for the orderings, and
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the representation.
