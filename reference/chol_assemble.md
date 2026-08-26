# The Cholesky Factor Behind a Free Vector

Assembles \\L\\ from the free vector: the diagonal is the exponential of
the first `dimension` values and the rest are placed below it, at the
positions
[`chol_positions()`](https://statmodels7.github.io/parameters7/reference/chol_positions.md)
recorded. Both
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
start here, and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
returns it unchanged.

## Usage

``` r
chol_assemble(s, eta)
```

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object, whose `dimension` and `param_params$positions` are read.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A `s@dimension` by `s@dimension` lower triangular numeric matrix with a
strictly positive diagonal, and no dimnames.

## Details

The exponential is applied by subsetting on `on_diagonal`, never through
[`ifelse()`](https://rdrr.io/r/base/ifelse.html), which evaluates both
branches over the whole vector. Here that would only exponentiate values
it then discards; the same shape in
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
would take a logarithm of below-diagonal entries that are free to be
negative, and warn about the `NaN`s it throws away.

## See also

[`chol_positions()`](https://statmodels7.github.io/parameters7/reference/chol_positions.md)
for the layout, and
[`param_factor.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_factor.LogCholeskyParam.md),
which returns this directly.
