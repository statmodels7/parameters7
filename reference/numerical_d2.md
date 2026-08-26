# Numerical Second Derivatives of a Parameter's Matrix

Estimates the \\d(d+1)/2\\ distinct second derivatives \\\partial^2 V /
\partial \eta_k \partial \eta_l\\, taking exactly one difference for
each. Which quantity is differenced depends on what the family already
supplies: the analytic first derivative where there is one, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
itself where there is not. This is the method
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
dispatches to when a family has not written its own.

## Usage

``` r
numerical_d2(s, eta)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list of `choose(s@n_free + 1, 2)` estimates keyed as
`param_tuple_names(s)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## Three routes, one layer each

Writing \\V(\eta)\\ for
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and \\\partial_k V\\ for an analytic first derivative, the component
\\(k, l)\\ is

\$\$\partial\_{kl} V \approx \frac{\partial_k V(\eta + h e_l) -
\partial_k V(\eta - h e_l)}{2h}\$\$

whenever
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
is analytic. Where it is not, an off-diagonal pair takes the four-point
mixed stencil

\$\$\partial\_{kl} V \approx \frac{V(\eta + h_k e_k + h_l e_l) -
V(\eta + h_k e_k - h_l e_l) - V(\eta - h_k e_k + h_l e_l) + V(\eta - h_k
e_k - h_l e_l)}{4 h_k h_l},\$\$

and a diagonal pair the three-point second difference. All three carry a
truncation error of order \\h^2\\.

## Why none of them nests

The rule the toolkit follows is that differences are never composed **in
the same variable**: a difference of a difference multiplies the error
of the first stage into the second, and a fourth derivative built that
way is noise. Differencing an analytic first derivative is one layer,
since only one stage is numerical. The four-point stencil differences
two *different* components, and two differences along different
coordinates commute into a single product stencil, so it is one layer as
well. The diagonal case uses the second-difference stencil directly and
never sees a first difference at all.

The steps are the order-2 ones, \\\varepsilon^{1/4}\max(1,
\|\eta_k\|)\\, about \\1.2 \times 10^{-4}\\ near the origin, since it is
a second derivative being estimated whichever route is taken.

## See also

[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
the generic this serves,
[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
for the order below, and
[`mixed_stencil()`](https://statmodels7.github.io/parameters7/reference/mixed_stencil.md),
which the third and fourth orders use.

## Examples

``` r
# A scalar matrix: every order is the matrix again, so the error shows.
s <- scalar_matrix(2)
numerical_d2(s, 0.3)
#> $`log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
max(abs(numerical_d2(s, 0.3)[[1]] - param_value(s, 0.3)))
#> [1] 7.642331e-12

# Against a family that writes its own second derivatives.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
max(abs(unlist(numerical_d2(q, eta)) - unlist(param_d2(q, eta))))
#> [1] 1.37506e-10
```
