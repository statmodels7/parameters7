# The Number of Diagonal Entries Each Free Value Owns

Returns \\m_k\\, the number of diagonal entries the \\k\\-th free value
controls: 1 for every value of a
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
and \\p\\ for the single value of a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).
It is the multiplicity that appears in the log-determinant's
derivatives, \\\log\|M\|\\ counting a shared value once per entry.

## Usage

``` r
diag_multiplicity(s)
```

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object, whose `param_params$shared`, `dimension` and `n_free` are
  read.

## Value

An integer vector of length `s@n_free` summing to `s@dimension`.

## See also

[`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md),
the same information per entry, and
[`param_dlogdet.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DiagMatrixParam.md),
which multiplies by it.
