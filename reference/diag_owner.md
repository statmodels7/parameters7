# The Free Value Each Diagonal Entry Belongs To

Returns, for each of the \\p\\ diagonal entries, the index into the free
vector of the value that controls it: `1:p` for a
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
and `rep(1, p)` for a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
whose one value controls them all. The derivative methods use it to
place a link derivative in the right entries.

## Usage

``` r
diag_owner(s)
```

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object, whose `param_params$shared` and `dimension` are read.

## Value

An integer vector of length `s@dimension`, with values in `1:s@n_free`.

## See also

[`diag_multiplicity()`](https://statmodels7.github.io/parameters7/reference/diag_multiplicity.md),
the same information counted the other way, and
[`diag_entries()`](https://statmodels7.github.io/parameters7/reference/diag_entries.md).
