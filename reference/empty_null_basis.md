# No Null Space

The empty basis a full-rank family declares. Every family whose value is
positive definite at every free vector passes this to its
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
constructor. The validator's shape rule, `dimension` by
`dimension - rank`, asks for exactly this when `rank` is `dimension`.

## Usage

``` r
empty_null_basis(dimension)
```

## Arguments

- dimension:

  The side of the matrix, a single integer.

## Value

A `dimension` by 0 numeric matrix.
[`ncol()`](https://rdrr.io/r/base/nrow.html) is 0, so a consumer that
loops over the null directions does nothing without a special case.
