# The Free Vector Behind a Matrix

The inverse map: given a matrix in the family's set, returns the free
vector that produces it.

## Usage

``` r
param_free(s, m, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- m:

  A symmetric numeric matrix of the parameter's dimension.

- ...:

  Passed to methods.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

Exact or refused, never obtained by optimisation. A generic
optimisation-based inverse would return a plausible \\\eta\\ for a
matrix outside the set, which is a wrong answer wearing the shape of a
right one. A family whose map has no closed-form inverse refuses
instead, and the base class refuses on behalf of any parameter that does
not implement this.

## See also

[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)

## Examples

``` r
s <- log_cholesky(2)
eta <- c(0.3, -0.2, 0.5)
param_free(s, param_value(s, eta))
#> log_L1 log_L2   L2.1 
#>    0.3   -0.2    0.5 
```
