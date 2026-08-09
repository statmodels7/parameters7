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

## Details

Let the tuple's distinct components be \\c_1, \dots, c_r\\ with
multiplicities \\m_1, \dots, m_r\\ summing to the order, and let
\\(o^{(m)}\_i, w^{(m)}\_i)\\ be the offsets and weights of the central
stencil for an \\m\\-th derivative in one variable. The estimate is
their tensor product,

\$\$\partial\_{c_1}^{m_1} \cdots \partial\_{c_r}^{m_r} M \approx
\frac{1}{\prod\_{j} h_j^{m_j}} \sum\_{i_1, \dots, i_r}
\left(\prod\_{j=1}^{r} w^{(m_j)}\_{i_j}\right) M\Bigl(\eta +
\textstyle\sum\_{j=1}^{r} o^{(m_j)}\_{i_j} h_j e\_{c_j}\Bigr),\$\$

evaluated at \\\prod_j (m_j + 1)\\ points of the map itself. Applying
one such stencil rather than differencing a lower-order numerical
derivative keeps the rounding at a single layer.

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
