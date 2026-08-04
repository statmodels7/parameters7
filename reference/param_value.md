# The Matrix a Parameter Produces

Maps the free vector \\\eta\\ to the matrix it parametrises.

## Usage

``` r
param_value(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A symmetric numeric matrix with `s@dimension` rows and columns.

## Details

This is the only generic a parameter must implement. Everything else in
the package has a numerical method registered on the
[`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
class and is therefore available from this one alone.

The generic validates the free vector before dispatching, so every
method, including one written outside the package, refuses a vector of
the wrong length or one containing a non-finite value. The free scale
has no boundary to reach, so a non-finite entry is a defect in the
caller rather than a point of the domain.

## See also

[`param_free`](https://statmodels7.github.io/parameters7/reference/param_free.md),
[`param_d1`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_logdet`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)

## Examples

``` r
s <- log_cholesky(2)
round(param_value(s, c(0, 0, 0.5)), 4)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25
```
