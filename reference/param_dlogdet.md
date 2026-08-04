# Gradient of the Log-Determinant

Returns \\\partial \log\|M\| / \partial \eta_k\\ for every free value,
or the same derivative of the log pseudo-determinant when the family is
rank deficient.

## Usage

``` r
param_dlogdet(s, eta, ...)
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

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

The identity behind it is \\\partial_k \log\|M\| = \mathrm{tr}(M^{-1}
\partial_k M)\\, with the Moore-Penrose inverse in place of \\M^{-1}\\
in the rank-deficient case. A closed form is therefore never an
independent claim: it must agree with
[`param_d1`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
through that identity, and
[`check_parameter`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
compares the two routes.

## See also

[`param_logdet`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
[`param_d2logdet`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)

## Examples

``` r
# for a scaled precision the derivative is the rank, whatever the scale
s <- scaled_matrix(crossprod(diff(diag(6), differences = 2)))
c(param_dlogdet(s, -3), param_dlogdet(s, 5), rank = s@rank)
#> log_scale log_scale      rank 
#>         4         4         4 
```
