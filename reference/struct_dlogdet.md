# Gradient of the Log-Determinant

Returns \\\partial \log\|M\| / \partial \eta_k\\ for every free value,
or the same derivative of the log pseudo-determinant when the family is
rank deficient.

## Usage

``` r
struct_dlogdet(s, eta, ...)
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

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

The identity behind it is \\\partial_k \log\|M\| = \mathrm{tr}(M^{-1}
\partial_k M)\\, with the Moore-Penrose inverse in place of \\M^{-1}\\
in the rank-deficient case. A closed form is therefore never an
independent claim: it must agree with
[`struct_dmatrix`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.md)
through that identity, and
[`check_covstruct`](https://statmodels7.github.io/covstructs7/reference/check_covstruct.md)
compares the two routes.

## See also

[`struct_logdet`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.md),
[`struct_d2logdet`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.md)

## Examples

``` r
# for a scaled precision the derivative is the rank, whatever the scale
s <- scaled_struct(crossprod(diff(diag(6), differences = 2)))
c(struct_dlogdet(s, -3), struct_dlogdet(s, 5), rank = s@rank)
#> scale scale  rank 
#>     4     4     4 
```
