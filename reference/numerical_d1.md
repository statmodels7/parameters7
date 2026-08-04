# Numerical First Derivatives of a Parameter's Matrix

One central difference of
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component of \\\eta\\.

## Usage

``` r
numerical_d1(s, eta)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A named list of symmetric matrices.

## See also

[`param_d1`](https://statmodels7.github.io/parameters7/reference/param_d1.md)

## Examples

``` r
numerical_d1(scalar_matrix(2), 0.3)
#> $scale
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
```
