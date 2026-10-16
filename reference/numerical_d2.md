# Numerical Second Derivatives of a Parameter's Matrix

Estimates the \\d(d+1)/2\\ distinct second derivatives \\\partial^2 V /
\partial \eta_k \partial \eta_l\\, taking exactly one difference for
each, through
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md).
Which quantity is differenced depends on what the family already
supplies: the analytic first derivative where there is one, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
itself where there is not. This is the method
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
dispatches to when a family has not written its own.

## Usage

``` r
numerical_d2(s, eta, accuracy = 4L)
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

A list of `choose(s@n_free + 1, 2)` estimates keyed as
`param_tuple_names(s)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## The two routes

Where
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
is analytic, the component \\(k, l)\\ is one central difference of
\\\partial_k V\\ along \\\eta_l\\. Where it is not, a diagonal component
is the central second difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and an off-diagonal one the tensor product of two first differences, in
\\\eta_k\\ and \\\eta_l\\. Differences are never composed in the same
variable, and two differences along different coordinates form one
product stencil, so every route is one layer.

The step of a difference of order \\r\\ is \\\varepsilon^{1/(r +
a)}\max(1, \|\eta_k\|)\\, with \\a\\ the accuracy: \\r = 1\\ on the
analytic route and \\r = 2\\ on the value. The truncation error is of
order \\h^a\\.

## See also

[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
the generic that this serves,
[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
for the order below, and
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md),
which builds the stencils.

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
#> [1] 1.063594e-13

# Against a family that writes its own second derivatives.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
max(abs(unlist(numerical_d2(q, eta)) - unlist(param_d2(q, eta))))
#> [1] 1.417533e-12
```
