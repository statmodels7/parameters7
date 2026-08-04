# Fourth Derivatives of a Parameter's Value

Returns the distinct fourth derivatives, keyed as
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, 4)`,
each shaped like the value.

## Usage

``` r
param_d4(s, eta, ...)
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

A named list, keyed as `param_tuple_names(s, 4)`.

## See also

[`param_d3`](https://statmodels7.github.io/parameters7/reference/param_d3.md)

## Examples

``` r
param_d4(scalar_matrix(2), 0.3)
#> $`log_scale:log_scale:log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
```
