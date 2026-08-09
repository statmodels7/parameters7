# Numerical Fourth Derivatives of a Parameter's Value

The order-four analogue of
[`numerical_d3`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md):
one stencil per component, applied directly to the map. Rounding is
amplified by the fourth power of the step, so this is accurate to
roughly four significant digits – a starting point, and the reason every
shipped family carries closed forms instead.

## Usage

``` r
numerical_d4(s, eta)
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
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)`(s, 4)`.

## Details

The tensor-product stencil is the one written out under
[`numerical_d3`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md),
with the multiplicities summing to four. Truncation is of order
\\h^{2}\\ and rounding of order \\\varepsilon / h^{4}\\, so the
attainable accuracy is roughly \\\varepsilon^{1/3}\\.

## See also

[`param_d4`](https://statmodels7.github.io/parameters7/reference/param_d4.md)

## Examples

``` r
numerical_d4(scalar_matrix(2), 0.3)
#> $`log_scale:log_scale:log_scale:log_scale`
#>          v1       v2
#> v1 1.349858 0.000000
#> v2 0.000000 1.349858
#> 
```
