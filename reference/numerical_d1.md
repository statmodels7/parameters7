# Numerical First Derivatives of a Parameter's Matrix

Estimates \\\partial V / \partial \eta_k\\ for every free value by one
three-point central difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component,

\$\$\partial_k V \approx \frac{V(\eta + h e_k) - V(\eta - h
e_k)}{2h},\$\$

at the step that
[`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
gives for a first derivative. This is the method
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
dispatches to when a family has not written its own, and it is exported
so that a family being developed can be compared against it.

## Usage

``` r
numerical_d1(s, eta)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list of `s@n_free` estimates named by `s@free_names`, each shaped like
the result of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and symmetrized for a matrix family.

## Details

The step is \\\varepsilon^{1/3}\max(1, \|\eta_k\|)\\, about \\6.1 \times
10^{-6}\\ near the origin, giving a truncation error of order \\h^2\\.
The example below compares the result with a closed form.

For a matrix family each estimate is symmetrized as `(m + t(m)) / 2`, as
the higher orders are, so that an asymmetry of rounding size in
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
does not reach the derivative. It costs \\2d\\ evaluations of the map.

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
#> [1] 7.642331e-12

# Against a family that writes its own.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
max(abs(unlist(numerical_d1(q, eta)) - unlist(param_d1(q, eta))))
#> [1] 6.8753e-11
```
