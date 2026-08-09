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

## Details

They follow from the second derivative of
[`param_d2logdet`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
by the same two rules, applied again: the trace is linear, and

\$\$\partial_m M^{-1} = -M^{-1}(\partial_m M)M^{-1}.\$\$

Every term of the result is therefore a trace of an alternating product
\\M^{-1}(\partial\_{I_1}M)M^{-1}(\partial\_{I_2}M)\cdots\\, one factor
per block of a partition of the index set, with the sign and the
multiplicity the two rules produce. The expansion is not transcribed:
the package differentiates
[`param_d1`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
through
[`param_d4`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
directly, and
[`check_parameter`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
holds the result against a numerical differentiation of the order below,
which shares none of its arithmetic.

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
