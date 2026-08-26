# The Factor of a Correlation Parameter, and Its Derivatives

Returns \\\partial^S L\\ for a multiset \\S\\ of free-value indices, or
`NULL` where that derivative is identically zero, which is most of the
time. The empty multiset gives \\L\\ itself.

## Usage

``` r
corr_dfactor(s, tb, ks)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object, whose `param_params$row` and `param_params$col` are read.

- tb:

  The tables of
  [`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md),
  evaluated at the point.

- ks:

  A multiset of free-value indices, possibly empty. Positions in
  `1:s@n_free`.

## Value

A `s@dimension` by `s@dimension` lower triangular numeric matrix, or
`NULL` where the derivative is identically zero.

## Details

Two rules make almost everything vanish. An entry of \\L\\ depends only
on the angles of **its own row**, so a multiset spanning two rows gives
the zero matrix. Within a row, an entry is a product over a prefix of
the angles, so it gives zero unless every differentiated angle appears
among its factors: a derivative in \\\theta\_{ij}\\ touches only the
entries of row \\i\\ from column \\j\\ onwards.

Note what this does **not** say about \\R\\. The factor's derivative
across two rows is zero; the derivative of \\R = LL^\top\\ across those
rows is not, \\R\_{ij}\\ being the inner product of two rows. See
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
for the support that does hold.

## See also

[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
for the input,
[`corr_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_derivative.md)
for the Leibniz sum that consumes it, and
[`chol_dfactor()`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md),
the same idea for the log-Cholesky family.
