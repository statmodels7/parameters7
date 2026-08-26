# The Symmetric Matrix Behind a Free Vector

Fills the lower triangle of \\S\\ with the free values, in the order
`param_params$positions` records, and mirrors it to the upper triangle.
No transformation is applied: the free values **are** the entries of
\\S\\, which is where this chart parts company with
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)'s.

## Usage

``` r
mlog_s(s, eta)
```

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object, whose `dimension` and `param_params$positions` are read.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A symmetric `s@dimension` by `s@dimension` numeric matrix, with no
dimnames and no constraint on its entries.

## See also

[`mlog_basis()`](https://statmodels7.github.io/parameters7/reference/mlog_basis.md)
for its derivative in one free value, and
[`param_value.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_value.MatrixLogParam.md),
which exponentiates it.
