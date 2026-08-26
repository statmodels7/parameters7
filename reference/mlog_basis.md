# The Basis Direction of One Free Value

Returns \\\partial S / \partial \eta_k\\, which is a constant matrix,
\\S\\ being linear in the free vector: a single 1 on the diagonal for a
diagonal free value, and a symmetric pair of 1s below and above the
diagonal for the rest. These are the directions the Frechet derivatives
contract against.

Because they do not depend on \\\eta\\, the whole nonlinearity of the
family sits in the exponential, and all four derivative orders are
contractions of the same fixed set of directions.

## Usage

``` r
mlog_basis(s, k)
```

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object, whose `dimension` and `param_params$positions` are read.

- k:

  The free-value index, a position in `1:s@n_free`.

## Value

A symmetric `s@dimension` by `s@dimension` numeric matrix with one or
two non-zero entries, all of them 1.

## See also

[`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md),
which rotates these by the eigenvectors of \\S\\, and
[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
which contracts them.
