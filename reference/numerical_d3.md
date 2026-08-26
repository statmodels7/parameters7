# Numerical Third Derivatives of a Parameter's Value

Estimates the distinct third derivatives by applying one product stencil
directly to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
per index tuple. A component the tuple repeats contributes the
one-dimensional stencil of the matching order, and a component appearing
once contributes a two-point central factor; the product is evaluated in
a single pass over the map. This is the method
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
dispatches to when a family has not written its own.

## Usage

``` r
numerical_d3(s, eta)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list of `choose(s@n_free + 2, 3)` estimates keyed as
`param_tuple_names(s, 3)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## The stencil

Let the tuple's distinct components be \\c_1, \dots, c_r\\ with
multiplicities \\m_1, \dots, m_r\\ summing to three, and let
\\(o^{(m)}\_i, w^{(m)}\_i)\\ be the offsets and weights of the central
stencil for an \\m\\-th derivative in one variable. The estimate is
their tensor product,

\$\$\partial\_{c_1}^{m_1} \cdots \partial\_{c_r}^{m_r} V \approx
\frac{1}{\prod\_{j} h_j^{m_j}} \sum\_{i_1, \dots, i_r}
\left(\prod\_{j=1}^{r} w^{(m_j)}\_{i_j}\right) V\Bigl(\eta +
\textstyle\sum\_{j=1}^{r} o^{(m_j)}\_{i_j} h_j e\_{c_j}\Bigr),\$\$

evaluated at one point per combination of nodes, with the zero-weight
nodes skipped. A tuple naming one component three times costs four
evaluations of the map, one naming a component twice and another once
six, and one naming three distinct components eight.

## Why it goes straight to the map

A lower-order numerical derivative is never differenced. The rounding of
a difference is amplified by \\h^{-k}\\, so two stages multiply their
errors and a third derivative built from a first and then a second is
far worse than one stencil of the order wanted. Going to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
directly keeps the rounding at a single layer whatever the family
implements, at the price of more evaluations of the map.

## Accuracy

Truncation is of order \\h^2\\ at a step of \\\varepsilon^{1/5}\max(1,
\|\eta_k\|)\\. Measured against the closed form of a \\2 \times 2\\
log-Cholesky covariance, the agreement is \\7 \times 10^{-6}\\ absolute
on entries of size 12, so about six digits.

## See also

[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md),
the generic this serves,
[`mixed_stencil()`](https://statmodels7.github.io/parameters7/reference/mixed_stencil.md),
which builds the product, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the order above.

## Examples

``` r
# A scalar matrix: every order is the matrix again, so the error shows.
s <- scalar_matrix(2)
numerical_d3(s, 0.3)
#> $`log_scale:log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
max(abs(numerical_d3(s, 0.3)[[1]] - param_value(s, 0.3)))
#> [1] 2.835123e-07

# Against a family that writes its own third derivatives.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
ana <- param_d3(q, eta)
c(gap = max(abs(unlist(numerical_d3(q, eta)) - unlist(ana))),
  scale = max(abs(unlist(ana))))
#>          gap        scale 
#> 7.352547e-06 1.193460e+01 
```
