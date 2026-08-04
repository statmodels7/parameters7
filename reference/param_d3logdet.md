# Third and Fourth Derivatives of the Log-Determinant

The higher derivatives of the log-(pseudo-)determinant, keyed as
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
of the matching order.

## Usage

``` r
param_d3logdet(s, eta, ...)

param_d4logdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A named numeric vector.

## See also

[`param_d2logdet`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)

## Examples

``` r
param_d3logdet(log_cholesky(2), c(0.1, -0.2, 0.4))
#> log_L1:log_L1:log_L1 log_L1:log_L1:log_L2   log_L1:log_L1:L2.1 
#>                    0                    0                    0 
#> log_L1:log_L2:log_L2   log_L1:log_L2:L2.1     log_L1:L2.1:L2.1 
#>                    0                    0                    0 
#> log_L2:log_L2:log_L2   log_L2:log_L2:L2.1     log_L2:L2.1:L2.1 
#>                    0                    0                    0 
#>       L2.1:L2.1:L2.1 
#>                    0 
```
