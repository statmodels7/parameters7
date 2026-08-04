# Second Derivatives of a Parameter's Value

Returns the \\d(d+1)/2\\ distinct second derivatives \\\partial^2 V /
\partial \eta_k \partial \eta_l\\, each shaped like the value.

## Usage

``` r
param_d2(s, eta, ...)
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

A named list keyed as `param_tuple_names(s)`.

## Details

The components are keyed by
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md),
generated from the same enumeration that produces the names rather than
by taking a name apart. Recovering an index by splitting a name on its
separator is the obvious route and it is wrong, because a free value
whose own label contains the separator splits into the wrong number of
pieces.

## See also

[`param_d1`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)

## Examples

``` r
s <- scaled_matrix(diag(2))
param_d2(s, 0)
#> $`log_scale:log_scale`
#>    v1 v2
#> v1  1  0
#> v2  0  1
#> 
```
