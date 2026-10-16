# Numerical First Derivatives of a Parameter's Matrix

Estimates \\\partial V / \partial \eta_k\\ for every free value by one
central difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component, through
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md).
This is the method
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
dispatches to when a family has not written its own, and it is exported
so that a family being developed can be compared against it.

## Usage

``` r
numerical_d1(s, eta, accuracy = 4L)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- accuracy:

  The order of accuracy of the central rule, a positive even whole
  number.

## Value

A list of `s@n_free` estimates named by `s@free_names`, each shaped like
the result of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and symmetrized for a matrix family.

## Details

At the default accuracy 4 the difference is the five-point rule
\\(V(\eta - 2h e_k) - 8V(\eta - h e_k) + 8V(\eta + h e_k) - V(\eta + 2h
e_k)) / 12h\\, at the step \\\varepsilon^{1/5}\max(1, \|\eta_k\|)\\,
about \\7.4 \times 10^{-4}\\ near the origin, with a truncation error of
order \\h^4\\. It costs \\4d\\ evaluations of the map; `accuracy = 2`
gives the three-point rule and \\2d\\.

For a matrix family each estimate is symmetrized as `(m + t(m)) / 2`, as
the higher orders are, so that an asymmetry of rounding size in
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
does not reach the derivative.

## See also

[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
the generic that this serves,
[`numerical_d2()`](https://statmodels7.github.io/parameters7/reference/numerical_d2.md),
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the higher orders, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
to find out whether a given family reaches this.

## Examples

``` r
# A scalar matrix is exp(eta) times the identity, so the derivative is the
# matrix itself and the error of the difference can be read off.
s <- scalar_matrix(2)
numerical_d1(s, 0.3)
#> $log_scale
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
max(abs(numerical_d1(s, 0.3)[[1]] - param_value(s, 0.3)))
#> [1] 1.063594e-13

# Against a family that writes its own.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
max(abs(unlist(numerical_d1(q, eta)) - unlist(param_d1(q, eta))))
#> [1] 7.087664e-13
```
