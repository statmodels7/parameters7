# Numerical Third Derivatives of a Parameter's Value

Estimates the distinct third derivatives, one tensor-product stencil per
index tuple, through
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md).
The stencil is applied to the highest analytic derivative below the
third that the family supplies, as
[`param_analytic()`](https://statmodels7.github.io/parameters7/reference/param_analytic.md)
reads it, and to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
where there is none. This is the method
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
dispatches to when a family has not written its own.

## Usage

``` r
numerical_d3(s, eta, accuracy = 4L)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- accuracy:

  The order of accuracy of the central rules, a positive even whole
  number.

## Value

A list of `choose(s@n_free + 2, 3)` estimates keyed as
`param_tuple_names(s, 3)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## The stencil

Let the tuple be \\\alpha\\ and the analytic quantity differenced be the
derivative \\\beta \le \alpha\\, taken from the components that the
tuple repeats most. The residual \\\alpha - \beta\\ has multiplicities
\\m_1, \dots, m_r\\ in the components \\c_1, \dots, c_r\\, and with
\\(o^{(m)}\_i, w^{(m)}\_i)\\ the offsets and weights of the central
stencil for an \\m\\-th derivative in one variable the estimate is their
tensor product,

\$\$\partial^\alpha V \approx \frac{1}{\prod\_{j} h_j^{m_j}} \sum\_{i_1,
\dots, i_r} \left(\prod\_{j=1}^{r} w^{(m_j)}\_{i_j}\right)
\partial^\beta V\Bigl(\eta + \textstyle\sum\_{j=1}^{r} o^{(m_j)}\_{i_j}
h_j e\_{c_j}\Bigr),\$\$

with the zero-weight nodes skipped and every point shared by several
tuples evaluated once. With an analytic
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
each tuple is one first difference of a second derivative.

## The stencil is never applied to a difference

A lower-order numerical derivative is never differenced. The rounding of
a difference is amplified by \\h^{-k}\\, so two stages multiply their
errors. The quantity differenced is either the map or an analytic
derivative, so the rounding stays at a single layer whatever the family
implements.

## Accuracy

Truncation is of order \\h^a\\ at a step of \\\varepsilon^{1/(r +
a)}\max(1, \|\eta_k\|)\\, with \\r\\ the order of the difference and
\\a\\ the accuracy, and rounding is of order \\\varepsilon / h^r\\.

## See also

[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md),
the generic that this serves,
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md),
which builds the stencils, and
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
#> [1] 1.063594e-13

# Against a family that writes its own third derivatives.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
ana <- param_d3(q, eta)
c(gap = max(abs(unlist(numerical_d3(q, eta)) - unlist(ana))),
  scale = max(abs(unlist(ana))))
#>          gap        scale 
#> 2.835066e-12 1.193460e+01 
```
