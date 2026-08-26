# The Diagonal Entries Behind a Free Vector

Applies the parameter's inverse link to the free vector and returns the
\\p\\ diagonal entries, recycling the single value across the whole
diagonal for a
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).
Every method of the family that needs the entries starts here, so the
shared and unshared cases branch in one place.

## Usage

``` r
diag_entries(s, eta)
```

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object, whose `param_params$link`, `param_params$shared` and
  `dimension` are read.

- eta:

  A numeric vector of free values, of length `s@n_free`: \\p\\ values
  for
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
  one for
  [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).

## Value

A numeric vector of length `s@dimension`, strictly positive under any
link
[`check_positive_link()`](https://statmodels7.github.io/parameters7/reference/check_positive_link.md)
admits.

## See also

[`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md)
for which free value each entry belongs to, and
[`diag_multiplicity()`](https://statmodels7.github.io/parameters7/reference/diag_multiplicity.md)
for how many entries each free value owns.
