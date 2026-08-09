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

Writing \\M(\eta)\\ for
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and \\\partial_k M\\ for its analytic first derivative, the entry \\(k,
l)\\ is

\$\$\partial\_{kl} M \approx \frac{\partial_k M(\eta + h e_l) -
\partial_k M(\eta - h e_l)}{2h},\$\$

one central difference of an analytic quantity. Where the first
derivatives are themselves numerical the four-point mixed stencil

\$\$\partial\_{kl} M \approx \frac{M(\eta + h e_k + h e_l) - M(\eta + h
e_k - h e_l) - M(\eta - h e_k + h e_l) + M(\eta - h e_k - h
e_l)}{4h^{2}}, \qquad k \ne l,\$\$

is used instead, with the three-point second difference on the diagonal.
Both are a single layer, and both carry a truncation error of order
\\h^{2}\\.

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
