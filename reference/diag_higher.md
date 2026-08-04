# Third and Fourth Derivatives of a Diagonal Matrix

The derivative components of orders three and four, each a diagonal
matrix with a single non-zero entry.

## Usage

``` r
diag_higher(s, eta, order)
```

## Arguments

- s:

  A
  [`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 3 or 4.

## Value

A named list of matrices, keyed by
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, order)`.

## Details

A diagonal family is separable, so a component is zero unless every
index of the tuple names the same free value, and it then carries the
corresponding derivative of the inverse link in the entry that value
owns. A shared free value owns every entry, which is what
[`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md)
answers.

## See also

[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
