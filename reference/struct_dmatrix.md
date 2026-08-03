# First Derivatives of a Structure's Matrix

Returns \\\partial M / \partial \eta_k\\ for every free value, as a list
of matrices.

## Usage

``` r
struct_dmatrix(s, eta, ...)
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

A list of `s@n_free` symmetric matrices, named by `s@free_names`.

## Details

A structure that registers no method for this generic gets the numerical
one of the
[`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
class, which applies a single central difference to
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md)
in each component.

## See also

[`struct_d2matrix`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.md),
[`struct_is_numerical`](https://statmodels7.github.io/covstructs7/reference/struct_is_numerical.md)

## Examples

``` r
s <- scaled_struct(diag(2))
struct_dmatrix(s, 0)
#> $scale
#>    v1 v2
#> v1  1  0
#> v2  0  1
#> 
```
