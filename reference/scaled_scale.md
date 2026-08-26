# The Scale Behind a Free Vector, and Its Derivatives

Returns the inverse link's value and its first two derivatives at the
single free value, so the family's
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
all read one call. For a fixed parameter, built with `link = NULL`, it
returns \\h = 1\\ and both derivatives 0, which makes those three
methods work unchanged: the matrix is \\P\\ and its derivatives are
empty lists.

## Usage

``` r
scaled_scale(s, eta)
```

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object, whose `param_params$link` is read.

- eta:

  A numeric vector of free values: one value, or `numeric(0)` for a
  fixed parameter.

## Value

A list with three single numbers, `h`, `d1` and `d2`.

## See also

[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the higher orders of \\\log h\\, and
[`param_value.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.ScaledMatrixParam.md),
the first caller.
