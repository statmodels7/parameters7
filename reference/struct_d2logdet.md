# Hessian of the Log-Determinant

Returns the distinct second derivatives of the log-determinant, or of
the log pseudo-determinant, keyed as
[`struct_pair_names`](https://statmodels7.github.io/covstructs7/reference/struct_pair_names.md).

## Usage

``` r
struct_d2logdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A named numeric vector, keyed as `struct_pair_names(s)`.

## See also

[`struct_dlogdet`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.md)

## Examples

``` r
struct_d2logdet(log_cholesky(2), c(0.2, -0.1, 0.4))
#> log_L1:log_L1 log_L2:log_L2     L2.1:L2.1 log_L1:log_L2   log_L1:L2.1 
#>             0             0             0             0             0 
#>   log_L2:L2.1 
#>             0 
```
