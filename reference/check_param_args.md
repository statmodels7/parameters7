# Validate the Arguments Shared by Every Matrix Constructor

Checks the two arguments every matrix family's constructor takes, and
returns the matrix side coerced to integer so the caller can store it in
the class's integer property. Called by
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
and
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md).

## Usage

``` r
check_param_args(dimension, role)
```

## Arguments

- dimension:

  The side of the matrix. Must be a single finite number, at least 1,
  equal to its own [`round()`](https://rdrr.io/r/base/Round.html). `0`,
  `2.5`, `c(1, 2)`, `"3"`, `Inf` and `NA` all throw
  `'dimension' must be a single positive integer.`

- role:

  The role label. Must be a single string, one of `"covariance"`,
  `"precision"` or `"either"`; anything else throws. The shipped
  constructors run
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) first, so their
  callers see [`match.arg()`](https://rdrr.io/r/base/match.arg.html)'s
  message and never this one. The check is here for a constructor
  written outside the package.

## Value

`dimension`, as a single integer.
