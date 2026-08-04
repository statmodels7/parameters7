# Third Derivatives of a Parameter's Value

Returns the distinct third derivatives \\\partial^3 V / \partial \eta_k
\partial \eta_l \partial \eta_m\\, keyed as
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, 3)`,
each shaped like the value.

## Usage

``` r
param_d3(s, eta, ...)
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

A named list, keyed as `param_tuple_names(s, 3)`.

## See also

[`param_d2`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d4`](https://statmodels7.github.io/parameters7/reference/param_d4.md)

## Examples

``` r
param_d3(scalar_matrix(2), 0.3)
#> $`log_scale:log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
```
