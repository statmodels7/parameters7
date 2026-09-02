# Validate the Argument Shared by Every Matrix Constructor

Checks the matrix side every matrix family's constructor takes, and
returns it coerced to integer so the caller can store it in the class's
integer property. Called by
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
check_param_args(dimension)
```

## Arguments

- dimension:

  The side of the matrix. Must be a single finite number, at least 1,
  equal to its own [`round()`](https://rdrr.io/r/base/Round.html). `0`,
  `2.5`, `c(1, 2)`, `"3"`, `Inf` and `NA` all throw
  `'dimension' must be a single positive integer.`

## Value

`dimension`, as a single integer.
