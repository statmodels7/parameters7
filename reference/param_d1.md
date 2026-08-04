# First Derivatives of a Parameter's Value

Returns \\\partial V / \partial \eta_k\\ for every free value, as a list
of objects each shaped like the value.

## Usage

``` r
param_d1(s, eta, ...)
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

A list of `s@n_free` symmetric matrices, named by `s@free_names`.

## Details

A parameter that registers no method for this generic gets the numerical
one of the
[`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
class, which applies a single central difference to
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component.

## See also

[`param_d2`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_is_numerical`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)

## Examples

``` r
s <- scaled_matrix(diag(2))
param_d1(s, 0)
#> $log_scale
#>    v1 v2
#> v1  1  0
#> v2  0  1
#> 
```
