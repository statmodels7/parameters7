# The Factor of a Log-Cholesky Parameter, and Its Derivatives

Returns \\\partial^S L\\ for a multiset \\S\\ of free-value indices, or
`NULL` where that derivative is identically zero. Almost all of them
are, and that is what keeps the Leibniz sum of
[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md)
short.

## Usage

``` r
chol_dfactor(s, l, ks)
```

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object, whose `dimension` and `param_params$positions` are read.

- l:

  The factor at the point, from
  [`chol_assemble()`](https://statmodels7.github.io/parameters7/reference/chol_assemble.md).

- ks:

  A multiset of free-value indices, possibly empty. Positions in
  `1:s@n_free`; the empty multiset gives \\L\\ itself.

## Value

A `s@dimension` by `s@dimension` numeric matrix with at most one
non-zero entry, or `l` unchanged for an empty `ks`, or `NULL` where the
derivative is identically zero.

## Details

Three cases exhaust it. An empty \\S\\ gives \\L\\ itself. A single
index gives a **single-entry matrix**: \\L\_{ii} E\_{ii}\\ for a
diagonal free value, whose parametrization is a logarithm, and
\\E\_{ij}\\ for one below the diagonal, which enters \\L\\ linearly. A
repeated index gives \\L\_{ii} E\_{ii}\\ again where every repetition
names the same diagonal value, every derivative of \\\exp(\eta_i)\\
being itself, and `NULL` otherwise.

So the only surviving higher derivatives are the pure diagonal ones, and
the single-entry shape is what the compiled kernel exploits: each
Leibniz term becomes one row, one column or one cell of a product, never
a matrix multiplication.

## See also

[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
the caller, and
[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md),
the generic Leibniz sum the R twin uses.
