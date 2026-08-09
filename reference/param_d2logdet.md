# Hessian of the Log-Determinant

Returns the distinct second derivatives of the log-determinant, or of
the log pseudo-determinant, keyed as
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md).

## Usage

``` r
param_d2logdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A named numeric vector, keyed as `param_tuple_names(s)`.

## Details

Differentiating \\\partial_k \log\|M\| = \mathrm{tr}(M^{-1}\partial_k
M)\\ once more, and using \\\partial_l M^{-1} = -M^{-1}(\partial_l
M)M^{-1}\\,

\$\$\partial\_{kl} \log\|M\| =
\mathrm{tr}\\\left(M^{-1}\partial\_{kl}M\right) -
\mathrm{tr}\\\left(M^{-1}(\partial_k M) M^{-1} (\partial_l
M)\right),\$\$

with the Moore-Penrose inverse in place of \\M^{-1}\\ when the family is
rank deficient. The second term is what makes the Hessian of the
log-determinant differ from the trace of the second derivative of the
matrix, and it is where a closed form is usually got wrong;
[`check_parameter`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
compares this route against
[`param_d2`](https://statmodels7.github.io/parameters7/reference/param_d2.md).

## See also

[`param_dlogdet`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)

## Examples

``` r
param_d2logdet(log_cholesky(2), c(0.2, -0.1, 0.4))
#> log_L1:log_L1 log_L2:log_L2     L2.1:L2.1 log_L1:log_L2   log_L1:L2.1 
#>             0             0             0             0             0 
#>   log_L2:L2.1 
#>             0 
```
