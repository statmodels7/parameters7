# Numerical Third Derivatives of a Parameter's Value

One product stencil applied directly to
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)
for each distinct index tuple: a single stencil on the map itself, never
a stencil on a lower-order numerical derivative, so the no-nesting rule
holds whatever the family implements. A repeated component uses the
matching higher-order one-dimensional factor; distinct components each
contribute a central two-point factor, and the product is evaluated in
one pass.

## Usage

``` r
numerical_d3(s, eta)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A named list keyed as
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, 3)`.

## See also

[`param_d3`](https://statmodels7.github.io/parameters7/reference/param_d3.md)

## Examples

``` r
numerical_d3(scalar_matrix(2), 0.3)
#> $`log_scale:log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
```
