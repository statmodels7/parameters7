# Second Derivatives of a Structure's Matrix

Returns the \\d(d+1)/2\\ distinct second derivatives \\\partial^2 M /
\partial \eta_k \partial \eta_l\\.

## Usage

``` r
struct_d2matrix(s, eta, ...)
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

A named list of symmetric matrices, keyed as `struct_pair_names(s)`.

## Details

The components are keyed by
[`struct_pair_names`](https://statmodels7.github.io/covstructs7/reference/struct_pair_names.md),
generated from the same enumeration that produces the names rather than
by taking a name apart. Recovering an index by splitting a name on its
separator is the obvious route and it is wrong, because a free value
whose own label contains the separator splits into the wrong number of
pieces.

## See also

[`struct_dmatrix`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.md),
[`struct_pair_names`](https://statmodels7.github.io/covstructs7/reference/struct_pair_names.md)

## Examples

``` r
s <- scaled_struct(diag(2))
struct_d2matrix(s, 0)
#> $`scale:scale`
#>    v1 v2
#> v1  1  0
#> v2  0  1
#> 
```
