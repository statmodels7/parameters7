# Numerical Second Derivatives of a Parameter's Matrix

One central difference of the analytic first derivatives where a
parameter supplies them, and a single mixed stencil on
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)
where it does not.

## Usage

``` r
numerical_d2(s, eta)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A named list of symmetric matrices, keyed as `param_tuple_names(s)`.

## Details

Differentiating the analytic first derivative costs one
finite-difference layer rather than two, which is the rule the whole
toolkit follows: never compose differences in the same variable. When
the first derivatives are themselves numerical the two differences act
on different components for an off-diagonal pair, so they commute into a
four-point mixed stencil rather than compounding; the diagonal pairs use
the three-point second-difference stencil directly, which is again one
layer and not two.

## See also

[`param_d2`](https://statmodels7.github.io/parameters7/reference/param_d2.md)

## Examples

``` r
numerical_d2(scalar_matrix(2), 0.3)
#> $`log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
```
