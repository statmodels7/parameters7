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

## See also

[`param_d4`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
